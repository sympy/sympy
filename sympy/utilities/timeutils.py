"""Simple tools for timing functions' execution, when IPython is not available."""

from __future__ import annotations

import os
import timeit
import math
from typing import TypeVar, Callable, ParamSpec, Any

_scales: list[float] = [1e0, 1e3, 1e6, 1e9]
_units: list[str] = ['s', 'ms', '\N{GREEK SMALL LETTER MU}s', 'ns']

def timed(func: Callable[[], object], setup: str = "pass", limit: int | None = None) -> tuple[int, float, float, str]:
    """Adaptively measure execution time of a function."""
    timer = timeit.Timer(func, setup=setup)
    repeat, number = 3, 1
    for _ in range(1, 10):
        if timer.timeit(number) >= 0.2:
            break
        elif limit is not None and number >= limit:
            break
        else:
            number *= 10
    time_taken = min(timer.repeat(repeat, number)) / number
    if time_taken > 0.0:
        order = min(-int(math.floor(math.log10(time_taken)) // 3), 3)
    else:
        order = 3
    return (number, time_taken, time_taken * _scales[order], _units[order])

# Code for doing inline timings of recursive algorithms.
def __do_timings() -> set[str]:
    env_str = os.getenv('SYMPY_TIMINGS', '')
    res_list = [x.strip() for x in env_str.split(',')]
    return set(res_list)



_do_timings: set[str] = __do_timings()
_timestack: list[Any] | None = None

def _print_timestack(stack: list[Any], level: int = 1) -> None:
    print('-' * level, '%.2f %s%s' % (stack[2], stack[0], stack[3]))
    for s in stack[1]:
        _print_timestack(s, level + 1)

P = ParamSpec('P')
R = TypeVar('R')

def timethis(name: str) -> Callable[[Callable[P, R]], Callable[P, R]]:
    def decorator(func: Callable[P, R]) -> Callable[P, R]:
        if name not in _do_timings:
            return func

        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            from time import time
            global _timestack
            oldtimestack = _timestack
            _timestack = [getattr(func, '__name__', 'unknown'), [], 0.0, args]

            t1 = time()
            r = func(*args, **kwargs)
            t2 = time()

            _timestack[2] = t2 - t1
            if oldtimestack is not None:
                oldtimestack[1].append(_timestack)  # type: ignore[index, attr-defined]
                _timestack = oldtimestack
            else:
                _print_timestack(_timestack)
                _timestack = None
            return r

        return wrapper
    return decorator
