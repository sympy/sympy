"""Tests for simple tools for timing functions' execution. """

import signal
import time
from sympy.testing.pytest import raises, skip
from sympy.utilities.timeutils import timed, timeout, TimeoutError, TimeoutSignal, TimeoutThreading

def test_timed():
    result = timed(lambda: 1 + 1, limit=100000)
    assert result[0] == 100000 and result[3] == "ns", str(result)

    result = timed("1 + 1", limit=100000)
    assert result[0] == 100000 and result[3] == "ns"



def test_timeout():
    def will_timeout():
        with timeout(1):
            time.sleep(2)

    def will_not_timeout():
        with timeout(duration_seconds=2):
            time.sleep(1)

    raises(TimeoutError, will_timeout)
    will_not_timeout()


def _test_Timeout(Cls):
    def func_fast(a, b):
        return 42*a+b

    def func_slow(a, b):
        time.sleep(2)
        return 17*a+b

    # gen is intentionally created already here.
    gen = Cls.map(lambda f, *args: f(*args), [func_fast, func_slow], [2, 3], [4, 5], seconds_per_item=1)

    def will_timeout():
        with Cls(duration_seconds=1):
            time.sleep(2)

    def will_not_timeout():
        with Cls(2):
            time.sleep(1)

    raises(TimeoutError, will_timeout)
    will_not_timeout()

    # we deliberately test gen several seconds after its creation.
    result = list(gen)
    assert result[0] == 42*2 + 4
    assert isinstance(result[1], TimeoutError)
    assert len(result) == 2

def test_TimeoutThreading():
    _test_Timeout(TimeoutThreading)

def test_TimeoutSignal():
    if not hasattr(signal, 'SIGALRM'):
        skip("ALARM signal not present.")
    _test_Timeout(TimeoutSignal)
