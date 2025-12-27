from __future__ import annotations
from functools import wraps
from typing import Any, Callable, List, TypeVar, Protocol, cast

T = TypeVar("T")
T_co = TypeVar("T_co", covariant=True)


class RecurrenceMemoFunc(Protocol[T_co]):
    """Protocol representing a memoized recurrence function."""

    def __call__(self, n: int) -> T_co: ...
    def cache_length(self) -> int: ...
    def fetch_item(self, x: Any) -> Any: ...


__all__ = [
    "recurrence_memo",
    "assoc_recurrence_memo",
]


def recurrence_memo(
    initial: List[T]
) -> Callable[[Callable[[int, List[T]], T]], RecurrenceMemoFunc[T]]:
    """
    Memo decorator for sequences defined by recurrence.

    Parameters
    ==========
    initial : List[T]
        Initial elements of the sequence.

    Returns
    =======
    Callable
        A decorator that converts a function into a callable matching the
        :class:`~sympy.utilities.memoization.RecurrenceMemoFunc` protocol.

    Notes
    =====
    - `T` is a type variable representing the element type of the sequence.

    Examples
    ========
    >>> from sympy.utilities.memoization import recurrence_memo
    >>> @recurrence_memo([1])  # 0! = 1
    ... def factorial(n, prev):
    ...     return n * prev[-1]
    >>> factorial(4)
    24
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
        return cast(RecurrenceMemoFunc[T], g)

    return decorator


def assoc_recurrence_memo(
    base_seq: Callable[[int], T]
) -> Callable[[Callable[[int, int, List[List[T]]], T]], Callable[[int, int], T]]:
    """
    Memo decorator for associated sequences defined by recurrence starting from base.

    Parameters
    ==========
    base_seq : Callable[[int], T]
        Function to generate base sequence elements.

    Returns
    =======
    Callable
        A decorator that converts a function into a callable generating
        associated sequences with memoization.

    Notes
    =====
    - Only works for P[n,0] = base_seq(n) cases.
    - Only works when m <= n.

    Examples
    ========
    >>> from sympy.utilities.memoization import assoc_recurrence_memo
    >>> base_seq = lambda n: 1
    >>> @assoc_recurrence_memo(base_seq)
    ... def assoc(n, m, cache):
    ...     return cache[n-1][m-1] + cache[n-1][m]
    >>> assoc(4, 2)
    6
    """
    cache: List[List[T]] = []

    def decorator(f: Callable[[int, int, List[List[T]]], T]) -> Callable[[int, int], T]:
        @wraps(f)
        def g(n: int, m: int) -> T:
            L = len(cache)
            if n < L:
                return cache[n][m]

            for i in range(L, n + 1):
                # get base sequence
                F_i0 = base_seq(i)
                F_i_cache = [F_i0]
                cache.append(F_i_cache)
                # generate assoc sequence
                for j in range(1, i + 1):
                    F_ij = f(i, j, cache)
                    F_i_cache.append(F_ij)

            return cache[n][m]

        return g

    return decorator
