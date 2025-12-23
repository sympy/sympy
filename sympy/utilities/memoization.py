from __future__ import annotations
from functools import wraps
from typing import Callable, List, Any


def recurrence_memo(
    initial: List[Any],
) -> Callable[..., Any]:
    """
    Memo decorator for sequences defined by recurrence.

    Parameters
    ==========
    initial : list
        The initial values of the sequence.

    Returns
    =======
    function
        The decorated recurrence function.

    Examples
    ========
    >>> from sympy.utilities.memoization import recurrence_memo
    >>> @recurrence_memo([1])
    ... def factorial(n, prev):
    ...     return n * prev[-1]
    >>> factorial(4)
    24
    >>> factorial(3)
    6
    >>> factorial.cache_length()
    5
    >>> factorial.fetch_item(slice(2, 4))
    [2, 6]
    """
    cache: List[Any] = initial

    def decorator(f: Callable[[int, List[Any]], Any]) -> Callable[..., Any]:
        @wraps(f)
        def g(n: int) -> Any:
            L = len(cache)
            if n < L:
                return cache[n]
            for i in range(L, n + 1):
                cache.append(f(i, cache))
            return cache[-1]

        g.cache_length = lambda: len(cache)  # type: ignore[attr-defined]
        g.fetch_item = lambda x: cache[x]    # type: ignore[attr-defined]
        return g

    return decorator


def assoc_recurrence_memo(
    base_seq: Callable[[int], Any],
) -> Callable[..., Any]:
    """
    Memo decorator for associated sequences defined by recurrence starting from base.

    Parameters
    ==========
    base_seq : callable
        Callable to get base sequence elements.

    Returns
    =======
    function
        The decorated associated recurrence function.

    Notes
    =====
    - Works only for Pn0 = base_seq(0) cases
    - Works only for m <= n cases
    """
    cache: List[List[Any]] = []

    def decorator(f: Callable[[int, int, List[List[Any]]], Any]) -> Callable[[int, int], Any]:
        @wraps(f)
        def g(n: int, m: int) -> Any:
            L = len(cache)
            if n < L:
                return cache[n][m]

            for i in range(L, n + 1):
                F_i0 = base_seq(i)
                F_i_cache = [F_i0]
                cache.append(F_i_cache)

                for j in range(1, i + 1):
                    F_ij = f(i, j, cache)
                    F_i_cache.append(F_ij)

            return cache[n][m]

        return g

    return decorator
