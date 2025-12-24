from __future__ import annotations
from functools import wraps
from typing import Any, Callable, List, TypeVar, Protocol, cast

_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True)

class RecurrenceMemoFunc(Protocol[_T_co]):
    def __call__(self, n: int) -> _T_co: ...
    def cache_length(self) -> int: ...
    def fetch_item(self, x) -> Any: ...

__all__ = ["recurrence_memo", "assoc_recurrence_memo", "RecurrenceMemoFunc"]

def recurrence_memo(
        initial: List[_T]
        ) -> Callable[
            [Callable[[int, List[_T]], _T]],RecurrenceMemoFunc[_T]
            ]:
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
    cache: List[_T] = initial

    def decorator(f: Callable[[int, List[_T]], _T]) -> RecurrenceMemoFunc[_T]:
        @wraps(f)
        def g(n: int) -> _T:
            L = len(cache)
            if n < L:
                return cache[n]
            for i in range(L, n + 1):
                cache.append(f(i, cache))
            return cache[-1]

        g.cache_length = lambda: len(cache)  # type: ignore[attr-defined]
        g.fetch_item = lambda x: cache[x]    # type: ignore[attr-defined]
        return cast(RecurrenceMemoFunc[_T], g)

    return decorator


def assoc_recurrence_memo(
        base_seq: Callable[[int], _T]
        ) -> Callable[
            [Callable[[int, int, List[List[_T]]], _T]],Callable[[int, int], _T]
              ]:
    """
    Memo decorator for associated sequences defined by recurrence starting from base

    base_seq(n) -- callable to get base sequence elements

    XXX works only for Pn0 = base_seq(0) cases
    XXX works only for m <= n cases
    """
    cache: List[List[_T]] = []

    def decorator(f: Callable[[int, int, List[List[_T]]], _T]) -> Callable[[int, int], _T]:
        @wraps(f)
        def g(n: int, m: int) -> _T:
            L = len(cache)
            if n < L:
                return cache[n][m]

            for i in range(L, n + 1):
                # get base sequence
                F_i0 = base_seq(i)
                F_i_cache = [F_i0]
                cache.append(F_i_cache)
                # XXX only works for m <= n cases
                # generate assoc sequence
                for j in range(1, i + 1):
                    F_ij = f(i, j, cache)
                    F_i_cache.append(F_ij)

            return cache[n][m]

        return g

    return decorator
