"""Simple tools for timing functions' execution, when IPython is not available. """
from __future__ import annotations

import timeit
import math
from typing import Any, Callable, TYPE_CHECKING, TypeVar

if TYPE_CHECKING:
    try:
        from typing import ParamSpec
    except ImportError:
        from typing_extensions import ParamSpec

    _P = ParamSpec('_P')
    _R = TypeVar('_R')

_scales = [1e0, 1e3, 1e6, 1e9]
_units = ['s', 'ms', '\N{GREEK SMALL LETTER MU}s', 'ns']


def timed(
    func: Callable[[], Any], setup: str = "pass", limit: int | None = None
) -> tuple[int, float, float, str]:
    """Adaptively measure execution time of a function. """
    timer = timeit.Timer(func, setup=setup)
    repeat, number = 3, 1

    for i in range(1, 10):
        if timer.timeit(number) >= 0.2:
            break
        elif limit is not None and number >= limit:
            break
        else:
            number *= 10

    time = min(timer.repeat(repeat, number)) / number

    if time > 0.0:
        order = min(-int(math.floor(math.log10(time)) // 3), 3)
    else:
        order = 3

    return (number, time, time*_scales[order], _units[order])


# Code for doing inline timings of recursive algorithms.

def __do_timings() -> set[str]:
    import os
    res = os.getenv('SYMPY_TIMINGS', '')
    res_list = [x.strip() for x in res.split(',')]
    return set(res_list)

_do_timings = __do_timings()
_timestack: list[Any] | None = None


def _print_timestack(stack: list[Any], level: int = 1) -> None:
    print('-'*level, '%.2f %s%s' % (stack[2], stack[0], stack[3]))
    for s in stack[1]:
        _print_timestack(s, level + 1)


def timethis(name: str) -> Callable[[Callable[_P, _R]], Callable[_P, _R]]:
    def decorator(func: Callable[_P, _R]) -> Callable[_P, _R]:
        if name not in _do_timings:
            return func

        def wrapper(*args: _P.args, **kwargs: _P.kwargs) -> _R:
            from time import time
            global _timestack
            oldtimestack = _timestack
            # Fixed legacy Python 2 func_name to Python 3 __name__
            _timestack = [func.__name__, [], 0, args]
            t1 = time()
            r = func(*args, **kwargs)
            t2 = time()
            _timestack[2] = t2 - t1
            if oldtimestack is not None:
                oldtimestack[1].append(_timestack)
                _timestack = oldtimestack
            else:
                _print_timestack(_timestack)
                _timestack = None
            return r
        return wrapper
    return decorator
    