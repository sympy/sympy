"""py.test hacks to support XFAIL/XPASS"""

import sys
import functools

try:
    import py
    from py.test import skip, raises
    USE_PYTEST = getattr(sys, '_running_pytest', False)
except ImportError:
    USE_PYTEST = False

if not USE_PYTEST:
    def raises(ExpectedException, code):
        """
        Tests that ``code`` raises the exception ``ExpectedException``.

        Does nothing if the right exception is raised, otherwise raises an
        AssertionError.

        Examples
        ========

        >>> from sympy.utilities.pytest import raises
        >>> raises(ZeroDivisionError, "1/0")
        >>> raises(ZeroDivisionError, "1/2")
        Traceback (most recent call last):
        ...
        AssertionError: DID NOT RAISE

        """
        if not isinstance(code, str):
            raise TypeError('raises() expects a code string for the 2nd argument.')
        frame = sys._getframe(1)
        loc = frame.f_locals.copy()
        try:
            exec code in frame.f_globals, loc
        except ExpectedException:
            return
        raise AssertionError("DID NOT RAISE")

    class XFail(Exception):
        pass

    class XPass(Exception):
        pass

    class Skipped(Exception):
        pass

    def XFAIL(func):
        def wrapper():
            try:
                func()
            except Exception, e:
                if sys.version_info[:2] < (2, 6):
                    message = getattr(e, 'message', '')
                else:
                    message = str(e)
                if message != "Timeout":
                    raise XFail(func.func_name)
                else:
                    raise Skipped("Timeout")
            raise XPass(func.func_name)

        wrapper = functools.update_wrapper(wrapper, func)
        return wrapper

    def skip(str):
        raise Skipped(str)


    def SKIP(reason):
        """Similar to :func:`skip`, but this is a decorator. """
        def wrapper(func):
            def func_wrapper():
                raise Skipped(reason)

            func_wrapper = functools.update_wrapper(func_wrapper, func)
            return func_wrapper

        return wrapper

    def slow(func):
        func._slow = True
        def func_wrapper():
            func()

        func_wrapper = functools.update_wrapper(func_wrapper, func)
        return func_wrapper

else:
    XFAIL = py.test.mark.xfail
    slow = py.test.mark.slow

    def SKIP(reason):
        def skipping(func):
            @functools.wraps(func)
            def inner(*args, **kwargs):
                skip(reason)
            return inner

        return skipping
