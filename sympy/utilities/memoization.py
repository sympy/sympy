
from __future__ import annotations

from functools import wraps
from typing import Callable, Protocol, TypeVar, cast, overload

T = TypeVar('T')


class _MemoizedSequence(Protocol[T]):
    def __call__(self, n: int) -> T: ...

    def cache_length(self) -> int: ...

    @overload
    def fetch_item(self, x: int) -> T: ...

    @overload
    def fetch_item(self, x: slice) -> list[T]: ...


def recurrence_memo(
    initial: list[T]
) -> Callable[[Callable[[int, list[T]], T]], _MemoizedSequence[T]]:

    cache: list[T] = initial

    def decorator(
        f: Callable[[int, list[T]], T]
    ) -> _MemoizedSequence[T]:

        @wraps(f)
        def g(n: int) -> T:
            L = len(cache)

            if n < L:
                return cache[n]

            for i in range(L, n + 1):
                cache.append(f(i, cache))

            return cache[-1]

        memoized_g = cast(_MemoizedSequence[T], g)

        memoized_g.cache_length = lambda: len(cache)
        memoized_g.fetch_item = lambda x: cache[x]

        return memoized_g

    return decorator


def assoc_recurrence_memo(
    base_seq: Callable[[int], T]
) -> Callable[
    [Callable[[int, int, list[list[T]]], T]],
    Callable[[int, int], T]
]:

    cache: list[list[T]] = []

    def decorator(
        f: Callable[[int, int, list[list[T]]], T]
    ) -> Callable[[int, int], T]:

        @wraps(f)
        def g(n: int, m: int) -> T:
            L = len(cache)

            if n < L:
                return cache[n][m]

            for i in range(L, n + 1):
                F_i0 = base_seq(i)
                F_i_cache: list[T] = [F_i0]
                cache.append(F_i_cache)

                for j in range(1, i + 1):
                    F_ij = f(i, j, cache)
                    F_i_cache.append(F_ij)

            return cache[n][m]

        return g

    return decorator
