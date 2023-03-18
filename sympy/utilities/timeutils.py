"""Simple tools for timing functions' execution, when IPython is not available. """

import contextlib
import os
import threading
import _thread
from typing import Callable
import queue
import time
import warnings
import signal
import timeit
import math


_scales = [1e0, 1e3, 1e6, 1e9]
_units = ['s', 'ms', '\N{GREEK SMALL LETTER MU}s', 'ns']


def timed(func, setup="pass", limit=None):
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

def __do_timings():
    import os
    res = os.getenv('SYMPY_TIMINGS', '')
    res = [x.strip() for x in res.split(',')]
    return set(res)

_do_timings = __do_timings()
_timestack = None


def _print_timestack(stack, level=1):
    print('-'*level, '%.2f %s%s' % (stack[2], stack[0], stack[3]))
    for s in stack[1]:
        _print_timestack(s, level + 1)


def timethis(name):
    def decorator(func):
        global _do_timings
        if name not in _do_timings:
            return func

        def wrapper(*args, **kwargs):
            from time import time
            global _timestack
            oldtimestack = _timestack
            _timestack = [func.func_name, [], 0, args]
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


class TimeoutError(Exception):
    pass


class _TimeoutContextManager(contextlib.AbstractContextManager):
    """Context manager to limit execution time.

    Examples
    --------
    >>> with TimeoutThreading(2) as t:
    ...     res = solve(expr, x)
    ...
    TimeoutError: Execution time exceeded 2 s.

    >>>
    >>> solutions = TimeoutThreading.map(solve, exprs, dep_vars, seconds_per_item=2)
    """

    def __init__(self, duration_seconds: int):
        self.duration_seconds = duration_seconds

    @classmethod
    def exception_type(cls) -> type[Exception]:
        return TimeoutError

    @classmethod
    def map(cls, func: Callable[..., ...], *args, seconds_per_item: int):
        for tup in zip(*args):
            try:
                with cls(duration_seconds=seconds_per_item):
                    yield func(*tup)
            except cls.exception_type() as e:
                yield e

    @property
    def elapsed(self) -> float:
        """Elapsed time in seconds."""
        return time.time() - self._t0

    @property
    def remaining(self) -> float:
        """Remaining time in seconds."""
        return self.duration_seconds - self.elapsed


    def exception_message(self) -> str:
        return f"Execution time exceeded {self.duration_seconds} s."

    def _construct_exception(self) -> Exception:
        return self.exception_type()(self.exception_message)

    def __enter__(self):
        self._t0 = time.time()


class TimeoutThreading(_TimeoutContextManager):
    def __enter__(self):
        super().__enter__()
        self._queue = queue.Queue()
        self.timer_thread = threading.Timer(self.duration_seconds, self._terminate)
        self.timer_thread.start()
        return self

    def _terminate(self):
        self._queue.put(self._construct_exception())
        _thread.interrupt_main()  # raises KeyboardInterrupt on the main thread.

    def __exit__(self, *exc):
        try:
            exc = self._queue.get(block=False)
        except queue.Empty:
            self.timer_thread.cancel()
            return False
        else:
            # there will be spurious KeyboardInterrupt exception in the traceback
            # if one really wanted to, it could be edited away, e.g.:
            # https://github.com/pallets/jinja/blob/9dad679695ccf84f06f8df897a6c5908e1963363/src/jinja2/debug.py#L14
            raise exc


class TimeoutSignal(_TimeoutContextManager):
    """Uses an ALARM signal (for POSIX compliant systems) to limit execution time."""
    def __enter__(self):
        super().__enter__()
        if os.name == 'nt':
            warnings.warn("Windows might not support ALARM signal, prefer TimeoutThreading instead.")
        signal.signal(signal.SIGALRM, self._terminate)
        signal.alarm(self.duration_seconds)
        return self

    def _terminate(self, signum, frame):
        assert signum == signal.SIGALRM
        raise self._construct_exception()

    def __exit__(self, *exc):
        signal.alarm(0)
        return False


default_timeout_context_manager_class = TimeoutThreading if os.name == 'nt' else TimeoutSignal


def timeout(duration_seconds: int):
    """Creates a context manager to limit the execution time of a code block.

    Parameters
    ----------
    duration_seconds: int
        Maximum execution time, after which a TimeoutError will be raised.

    Examples
    --------
    >>> with timeout(3) as t:
    ...     print("you will see this")
    ...     i = 0
    ...     for _ in range(10**20):
    ...         i += 1
    ...         if i % 10**7 == 0:
    ...             print(t.remaining)
    ...     print("but never this")  # DOCTEST: +ignore
    you will see this
    3.0
    2.21881103515625
    1.4843239784240723
    0.8278305530548096
    0.12401294708251953
    ---------------------------------------------------------------------------
    KeyboardInterrupt                         Traceback (most recent call last)
    ...
    During handling of the above exception, another exception occurred:

    TimeoutError                              Traceback (most recent call last)
    ...
    TimeoutError: Execution time exceeded 3 s.


    """
    return default_timeout_context_manager_class(duration_seconds)
