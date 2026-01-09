from functools import wraps
from typing import (
    TypeVar,
    Callable,
    List,
    Protocol,
    cast,
    Union,
    TYPE_CHECKING,
)



T = TypeVar("T")   # Elements of the main sequence
U = TypeVar("U")   # Elements of associated sequences



class RecurrenceFunc(Protocol[T]):
    def __call__(self, n: int) -> T: ...
    cache_length: Callable[[], int]
    fetch_item: Callable[[Union[int, slice]], Union[T, List[T]]]

if TYPE_CHECKING:
    recurrence_memo: Callable[
        [List[T]],
        Callable[[Callable[[int, List[T]], T]], RecurrenceFunc[T]],
    ]

    assoc_recurrence_memo: Callable[
        [Callable[[int], U]],
        Callable[
            [Callable[[int, int, List[List[U]]], U]],
            Callable[[int, int], U],
        ],
    ]

def recurrence_memo(initial):

    """
    Memo decorator for sequences defined by recurrence.

    Returns
    =======
    callable
        Function that can be called with an integer argument ``n``.

        The returned function has additional attributes:

        * ``cache_length()`` - return the current cache size
        * ``fetch_item(x)`` - return cached item(s) for an index or slice

    Examples
    ========

    >>> from sympy.utilities.memoization import recurrence_memo
    >>> @recurrence_memo([1])  # 0! = 1
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
    cache: List[T] = initial.copy()

    def decorator(f: Callable[[int, List[T]], T]) -> RecurrenceFunc[T]:
        @wraps(f)
        def g(n: int) -> T:
            L = len(cache)
            if n < L:
                return cache[n]
            for i in range(L, n + 1):
                cache.append(f(i, cache))
            return cache[-1]

        rf = cast(RecurrenceFunc[T], g)
        rf.cache_length = lambda: len(cache)
        rf.fetch_item = lambda x: cache[x]
        return rf

    return decorator

def assoc_recurrence_memo(base_seq):
    """
    Memo decorator for associated sequences defined by recurrence starting from base.

    base_seq(n) -- callable to get base sequence elements

    XXX works only for Pn0 = base_seq(0) cases
    XXX works only for m <= n cases
    """

    cache: List[List[U]] = []

    def decorator(f: Callable[[int, int, List[List[U]]], U]) -> Callable[[int, int], U]:
        @wraps(f)
        def g(n: int, m: int) -> U:
            L = len(cache)
            if n < L:
                return cache[n][m]

            for i in range(L, n + 1):
                F_i0 = base_seq(i)
                F_i_cache: List[U] = [F_i0]
                cache.append(F_i_cache)

                for j in range(1, i + 1):
                    F_ij = f(i, j, cache)
                    F_i_cache.append(F_ij)

            return cache[n][m]

        return g

    return decorator
