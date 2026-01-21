from functools import wraps
from typing import Callable, List, Union, cast, TYPE_CHECKING


if TYPE_CHECKING:
    from typing import Protocol, TypeVar

    T = TypeVar("T")

    class _RecurrenceFunc(Protocol[T]):
        def __call__(self, n: int) -> T: ...
        cache_length: Callable[[], int]
        fetch_item: Callable[[Union[int, slice]], Union[T, List[T]]]
else:
    # Runtime-safe fallbacks (avoid Protocol/TypeVar at runtime)
    T = object

    class _RecurrenceFunc:
        pass


def recurrence_memo(initial):
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
    cache = list(initial)

    def decorator(f):
        @wraps(f)
        def g(n):
            L = len(cache)
            if n < L:
                return cache[n]
            for i in range(L, n + 1):
                cache.append(f(i, cache))
            return cache[-1]

        rf = cast(_RecurrenceFunc, g)
        rf.cache_length = lambda: len(cache)
        rf.fetch_item = lambda x: cache[x]
        return rf

    return decorator


def assoc_recurrence_memo(base_seq):
    """
    Memo decorator for associated sequences defined by recurrence starting from base

    base_seq(n) -- callable to get base sequence elements

    XXX works only for Pn0 = base_seq(0) cases
    XXX works only for m <= n cases
    """
    cache = []

    def decorator(f):
        @wraps(f)
        def g(n, m):
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
