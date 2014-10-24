import time
from sympy.core.compatibility import TimeoutError
from functools import partial


def make_timeout_callback(timeout):
    """Creates a timeout callback function. Function will return True if
    the time is up, otherwise returns False."""

    TIME = time.time() + timeout
    def callback():
        """Return true if time is up, else False"""
        if time.time() > TIME:
            return True
        else:
            return False
    return callback

no_timeout = lambda: False

def init_timeout(arg):
    """Initiliaze the timeout callback.

    Parameters:
    -----------
    arg
        arg can be either a time in seconds, an existing callback function,
        or None. If None, returns a callback that always returns False.

    Returns:
    --------
    tcb : callable
        Timeout callback. Returns True if the time is up, else False.
    """

    if arg is None:
        return no_timeout
    elif callable(arg):
        return arg
    else:
        return make_timeout_callback(arg)

def check_timeout(tcb):
    """Throws TimeoutError if timeout callback is expired."""
    if tcb and tcb():
        raise TimeoutError()

def safe_partial_tcb(func, tcb):
    """Apply partial to a rule only if that rule has a tcb as a kwarg"""
    if 'tcb' in func.__code__.co_varnames:
        return partial(func, tcb=tcb)
    else:
        return func
