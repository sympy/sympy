from __future__ import annotations
from functools import wraps
from typing import Any, Callable, TypeVar

_V = TypeVar("_V")


def recurrence_memo(initial: list[_V]) -> Callable[[Callable[[int, list[_V]], _V]], Any]:
    """
    Memo decorator for sequences defined by recurrence

    Examples
    ========

    >>> from sympy.utilities.memoization import recurrence_memo
    >>> @recurrence_memo([1]) # 0! = 1
    ... def factorial(n, prev):
    ...     return n * prev[-1]
    >>> factorial(4)
    24
    >>> factorial(3) # use cache values
    6
    >>> factorial.cache_length() # cache length can be obtained
    5
    >>> factorial.fetch_item(slice(2, 4))
    [2, 6]

    """
    cache = initial

    def decorator(f: Callable[[int, list[_V]], _V]) -> Any:
        @wraps(f)
        def g(n: int) -> _V:
            L = len(cache)
            if n < L:
                return cache[n]
            for i in range(L, n + 1):
                cache.append(f(i, cache))
            return cache[-1]
        
        # We ignore type checks here because we're dynamically setting attributes on a function,
        # which standard 'Callable' doesn't easily support typing for without a Protocol.
        g.cache_length = lambda: len(cache)  # type: ignore[attr-defined]
        g.fetch_item = lambda x: cache[x]  # type: ignore[attr-defined]
        return g
    return decorator


def assoc_recurrence_memo(base_seq: Callable[[int], _V]) -> Callable[[Callable[[int, int, list[list[_V]]], _V]], Any]:
    """
    Memo decorator for associated sequences defined by recurrence starting from base

    base_seq(n) -- callable to get base sequence elements

    XXX works only for Pn0 = base_seq(0) cases
    XXX works only for m <= n cases
    """

    cache: list[list[_V]] = []

    def decorator(f: Callable[[int, int, list[list[_V]]], _V]) -> Any:
        @wraps(f)
        def g(n: int, m: int) -> _V:
            L = len(cache)
            if n < L:
                return cache[n][m]

            for i in range(L, n + 1):
                # get base sequence
                F_i0 = base_seq(i)
                F_i_cache: list[_V] = [F_i0]
                cache.append(F_i_cache)

                # XXX only works for m <= n cases
                # generate assoc sequence
                for j in range(1, i + 1):
                    F_ij = f(i, j, cache)
                    F_i_cache.append(F_ij)

            return cache[n][m]

        return g
    return decorator

