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
    def raises(expectedException, code=None):
        """
        Tests that ``code`` raises the exception ``expectedException``.

        ``code`` may be a callable, such as a lambda expression or function
        name.

        If ``code`` is not given or None, ``raises`` will return a context
        manager for use in ``with`` statements; the code to execute then
        comes from the scope of the ``with``. The calling module must have
        ``from __future__ import with_statement`` as its first code line
        for compatibility with Python 2.5 in that case.

        raises does nothing if the callable raises the expected exception,
        otherwise it raises an AssertionError.

        Examples
        ========

        >>> from sympy.utilities.pytest import raises

        >>> raises(ZeroDivisionError, lambda: 1/0)
        >>> raises(ZeroDivisionError, lambda: 1/2)
        Traceback (most recent call last):
        ...
        AssertionError: DID NOT RAISE

        (Python 2.5's doctest cannot support with statements due to a bug in
        its __import__ handling. That's why the following examples are
        excluded from doctesting; the doctest: annotations can go once SymPy
        stops Python 2.5 support.)

        >>> with raises(ZeroDivisionError): # doctest: +SKIP
        ...     n = 1/0
        >>> with raises(ZeroDivisionError): # doctest: +SKIP
        ...     n = 1/2
        Traceback (most recent call last):
        ...
        AssertionError: DID NOT RAISE

        Note that you cannot test multiple statements via
        ``with raises``:

        >>> with raises(ZeroDivisionError): # doctest: +SKIP
        ...     n = 1/0    # will execute and raise, aborting the ``with``
        ...     n = 9999/0 # never executed

        This is just what ``with`` is supposed to do: abort the
        contained statement sequence at the first exception and let
        the context manager deal with the exception.

        To test multiple statements, you'll need a separate ``with``
        for each:

        >>> with raises(ZeroDivisionError): # doctest: +SKIP
        ...     n = 1/0    # will execute and raise
        ... with raises(ZeroDivisionError):
        ...     n = 9999/0 # will also execute and raise

        """
        if code is None:
            return RaisesContext(expectedException)
        elif callable(code):
            try:
                code()
            except expectedException:
                return
            raise AssertionError("DID NOT RAISE")
        elif isinstance(code, str):
            raise TypeError(
                '\'raises(xxx, "code")\' has been phased out; '
                'change \'raises(xxx, "expression")\' '
                'to \'raises(xxx, lambda: expression)\', '
                '\'raises(xxx, "statement")\' '
                'to \'with raises(xxx): statement\'')
        else:
            raise TypeError(
                'raises() expects a callable for the 2nd argument.')

    class RaisesContext(object):
        def __init__(self, expectedException):
            self.expectedException = expectedException
        def __enter__(self):
            return None
        def __exit__(self, exc_type, exc_value, traceback):
            if exc_type is None:
                raise AssertionError("DID NOT RAISE")
            return issubclass(exc_type, self.expectedException)

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
