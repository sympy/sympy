from __future__ import annotations
from functools import wraps
from typing import Callable, List, TypeVar, Protocol, Any

T = TypeVar("T")
T_co = TypeVar("T_co", covariant=True)


class RecurrenceMemoFunc(Protocol[T_co]):
    def __call__(self, n: int) -> T_co: ...
    def cache_length(self) -> int: ...
    def fetch_item(self, x: Any): ...


def recurrence_memo(
    initial: List[T],
) -> Callable[[Callable[[int, List[T]], T]], RecurrenceMemoFunc[T]]:
    """
    Memo decorator for sequences defined by recurrence

    Examples
    ========

    >>> from sympy.utilities.memoization import recurrence_memo
    >>> @recurrence_memo([1])  # 0! = 1
    ... def factorial(n, prev):
    ...     return n * prev[-1]
    >>> factorial(4)
    24
    >>> factorial(3)  # use cache values
    6
    >>> factorial.cache_length()  # cache length can be obtained
    5
    >>> factorial.fetch_item(slice(2, 4))
    [2, 6]
    """
    cache: List[T] = initial

    def decorator(f: Callable[[int, List[T]], T]) -> RecurrenceMemoFunc[T]:
        @wraps(f)
        def g(n: int) -> T:
            L = len(cache)
            if n < L:
                return cache[n]
            for i in range(L, n + 1):
                cache.append(f(i, cache))
            return cache[-1]

        g.cache_length = lambda: len(cache)  # type: ignore[attr-defined]
        g.fetch_item = lambda x: cache[x]    # type: ignore[attr-defined]
        return g  # type: ignore[return-value]

    return decorator


def assoc_recurrence_memo(
    base_seq: Callable[[int], T],
) -> Callable[[Callable[[int, int, List[List[T]]], T]], Callable[[int, int], T]]:
    """
    Memo decorator for associated sequences defined by recurrence starting from base

    base_seq(n) -- callable to get base sequence elements

    XXX works only for Pn0 = base_seq(0) cases
    XXX works only for m <= n cases
    """
    cache: List[List[T]] = []

    def decorator(f: Callable[[int, int, List[List[T]]], T]) -> Callable[[int, int], T]:
        @wraps(f)
        def g(n: int, m: int) -> T:
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
