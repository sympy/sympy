"""py.test hacks to support XFAIL/XPASS"""

import sys
import functools
import os
import contextlib
import warnings
from typing import Any, Callable

from sympy.utilities.exceptions import SymPyDeprecationWarning

ON_TRAVIS = os.getenv('TRAVIS_BUILD_NUMBER', None)

try:
    import pytest
    USE_PYTEST = getattr(sys, '_running_pytest', False)
except ImportError:
    USE_PYTEST = False


raises: Callable[[Any, Any], Any]
XFAIL: Callable[[Any], Any]
skip: Callable[[Any], Any]
SKIP: Callable[[Any], Any]
slow: Callable[[Any], Any]
nocache_fail: Callable[[Any], Any]


if USE_PYTEST:
    raises = pytest.raises
    warns = pytest.warns
    skip = pytest.skip
    XFAIL = pytest.mark.xfail
    SKIP = pytest.mark.skip
    slow = pytest.mark.slow
    nocache_fail = pytest.mark.nocache_fail
    from _pytest.outcomes import Failed

else:
    # Not using pytest so define the things that would have been imported from
    # there.

    # _pytest._code.code.ExceptionInfo
    class ExceptionInfo:
        def __init__(self, value):
            self.value = value

        def __repr__(self):
            return "<ExceptionInfo {!r}>".format(self.value)


    def raises(expectedException, code=None):
        """
        Tests that ``code`` raises the exception ``expectedException``.

        ``code`` may be a callable, such as a lambda expression or function
        name.

        If ``code`` is not given or None, ``raises`` will return a context
        manager for use in ``with`` statements; the code to execute then
        comes from the scope of the ``with``.

        ``raises()`` does nothing if the callable raises the expected exception,
        otherwise it raises an AssertionError.

        Examples
        ========

        >>> from sympy.testing.pytest import raises

        >>> raises(ZeroDivisionError, lambda: 1/0)
        <ExceptionInfo ZeroDivisionError(...)>
        >>> raises(ZeroDivisionError, lambda: 1/2)
        Traceback (most recent call last):
        ...
        Failed: DID NOT RAISE

        >>> with raises(ZeroDivisionError):
        ...     n = 1/0
        >>> with raises(ZeroDivisionError):
        ...     n = 1/2
        Traceback (most recent call last):
        ...
        Failed: DID NOT RAISE

        Note that you cannot test multiple statements via
        ``with raises``:

        >>> with raises(ZeroDivisionError):
        ...     n = 1/0    # will execute and raise, aborting the ``with``
        ...     n = 9999/0 # never executed

        This is just what ``with`` is supposed to do: abort the
        contained statement sequence at the first exception and let
        the context manager deal with the exception.

        To test multiple statements, you'll need a separate ``with``
        for each:

        >>> with raises(ZeroDivisionError):
        ...     n = 1/0    # will execute and raise
        >>> with raises(ZeroDivisionError):
        ...     n = 9999/0 # will also execute and raise

        """
        if code is None:
            return RaisesContext(expectedException)
        elif callable(code):
            try:
                code()
            except expectedException as e:
                return ExceptionInfo(e)
            raise Failed("DID NOT RAISE")
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

    class RaisesContext:
        def __init__(self, expectedException):
            self.expectedException = expectedException

        def __enter__(self):
            return None

        def __exit__(self, exc_type, exc_value, traceback):
            if exc_type is None:
                raise Failed("DID NOT RAISE")
            return issubclass(exc_type, self.expectedException)

    class XFail(Exception):
        pass

    class XPass(Exception):
        pass

    class Skipped(Exception):
        pass

    class Failed(Exception):  # type: ignore
        pass

    def XFAIL(func):
        def wrapper():
            try:
                func()
            except Exception as e:
                message = str(e)
                if message != "Timeout":
                    raise XFail(func.__name__)
                else:
                    raise Skipped("Timeout")
            raise XPass(func.__name__)

        wrapper = functools.update_wrapper(wrapper, func)
        return wrapper

    def skip(str):
        raise Skipped(str)

    def SKIP(reason):
        """Similar to ``skip()``, but this is a decorator. """
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
        func_wrapper.__wrapped__ = func
        return func_wrapper

    def nocache_fail(func):
        "Dummy decorator for marking tests that fail when cache is disabled"
        return func

    @contextlib.contextmanager
    def warns(warningcls, *, match=''):
        '''Like raises but tests that warnings are emitted.

        >>> from sympy.testing.pytest import warns
        >>> import warnings

        >>> with warns(UserWarning):
        ...     warnings.warn('deprecated', UserWarning)

        >>> with warns(UserWarning):
        ...     pass
        Traceback (most recent call last):
        ...
        Failed: DID NOT WARN. No warnings of type UserWarning\
        was emitted. The list of emitted warnings is: [].
        '''
        # Absorbs all warnings in warnrec
        with warnings.catch_warnings(record=True) as warnrec:
            # Hide all warnings but make sure that our warning is emitted
            warnings.simplefilter("ignore")
            warnings.filterwarnings("always", match, warningcls)
            # Now run the test
            yield

        # Raise if expected warning not found
        if not any(issubclass(w.category, warningcls) for w in warnrec):
            msg = ('Failed: DID NOT WARN.'
                   ' No warnings of type %s was emitted.'
                   ' The list of emitted warnings is: %s.'
                   ) % (warningcls, [w.message for w in warnrec])
            raise Failed(msg)


def _both_exp_pow(func):
    """
    Decorator used to run the test twice: the first time `e^x` is represented
    as ``Pow(E, x)``, the second time as ``exp(x)`` (exponential object is not
    a power).

    This is a temporary trick helping to manage the elimination of the class
    ``exp`` in favor of a replacement by ``Pow(E, ...)``.
    """
    from sympy.core.parameters import _exp_is_pow

    def func_wrap():
        with _exp_is_pow(True):
            func()
        with _exp_is_pow(False):
            func()

    wrapper = functools.update_wrapper(func_wrap, func)
    return wrapper


@contextlib.contextmanager
def warns_deprecated_sympy():
    '''Shorthand for ``warns(SymPyDeprecationWarning)``

    This is the recommended way to test that ``SymPyDeprecationWarning`` is
    emitted for deprecated features in SymPy. To test for other warnings use
    ``warns``. To suppress warnings without asserting that they are emitted
    use ``ignore_warnings``.

    >>> from sympy.testing.pytest import warns_deprecated_sympy
    >>> from sympy.utilities.exceptions import SymPyDeprecationWarning
    >>> with warns_deprecated_sympy():
    ...     SymPyDeprecationWarning("Don't use", feature="old thing",
    ...         deprecated_since_version="1.0", issue=123).warn()

    >>> with warns_deprecated_sympy():
    ...     pass
    Traceback (most recent call last):
    ...
    Failed: DID NOT WARN. No warnings of type \
    SymPyDeprecationWarning was emitted. The list of emitted warnings is: [].
    '''
    with warns(SymPyDeprecationWarning):
        yield

@contextlib.contextmanager
def ignore_warnings(warningcls):
    '''Context manager to suppress warnings during tests.

    This function is useful for suppressing warnings during tests. The warns
    function should be used to assert that a warning is raised. The
    ignore_warnings function is useful in situation when the warning is not
    guaranteed to be raised (e.g. on importing a module) or if the warning
    comes from third-party code.

    When the warning is coming (reliably) from SymPy the warns function should
    be preferred to ignore_warnings.

    >>> from sympy.testing.pytest import ignore_warnings
    >>> import warnings

    Here's a warning:

    >>> with warnings.catch_warnings():  # reset warnings in doctest
    ...     warnings.simplefilter('error')
    ...     warnings.warn('deprecated', UserWarning)
    Traceback (most recent call last):
      ...
    UserWarning: deprecated

    Let's suppress it with ignore_warnings:

    >>> with warnings.catch_warnings():  # reset warnings in doctest
    ...     warnings.simplefilter('error')
    ...     with ignore_warnings(UserWarning):
    ...         warnings.warn('deprecated', UserWarning)

    (No warning emitted)
    '''
    # Absorbs all warnings in warnrec
    with warnings.catch_warnings(record=True) as warnrec:
        # Make sure our warning doesn't get filtered
        warnings.simplefilter("always", warningcls)
        # Now run the test
        yield

    # Reissue any warnings that we aren't testing for
    for w in warnrec:
        if not issubclass(w.category, warningcls):
            warnings.warn_explicit(w.message, w.category, w.filename, w.lineno)
