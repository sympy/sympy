from functools import wraps
import contextlib
import time

class _global_timeout(list):
    def __setitem__(self, key, value):
        super(_global_timeout, self).__setitem__(key, value)

TIMEOUT = _global_timeout([None])


@contextlib.contextmanager
def timeout_after(sec):
    """Context manager for controlling timeouts on sympy function calls.

    Raises a ``TimeoutError`` if the contained codeblock takes longer than
    `sec`. Note that calls not decorated with `timeoutable` will block.

    Nested `timeout_after` managers will not change the timeout. Once set at
    the top level, it will remain until the top level context is released.

    Parameters:
    -----------
    sec : number
        Maximum time in seconds that the containing code block can take.
    """
    nested = TIMEOUT[0] is not None
    try:
        if not nested and sec is not None:
            TIMEOUT[0] = time.time() + sec
        yield
    finally:
        if not nested:
            TIMEOUT[0] = None


def timeoutable(func):
    """Decorator indicating that the containing function allows timeouts to
    occur. To be used with `timeout_after` context-manager."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        if not TIMEOUT[0] or TIMEOUT[0] > time.time():
            return func(*args, **kwargs)
        else:
            raise TimeoutError()
    return wrapper
