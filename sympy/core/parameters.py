"""Thread-safe global parameters"""

from .cache import clear_cache
from contextlib import contextmanager
from threading import local

class _global_parameters(local):
    """
    Thread-local global parameters.

    Explanation
    ===========

    This class generates thread-local container for SymPy's global parameters.
    Every global parameters must be passed as keyword argument when generating
    its instance.
    A variable, `global_parameters` is provided as default instance for this class.

    WARNING! Although the global parameters are thread-local, SymPy's cache is not
    by now.
    This may lead to undesired result in multi-threading operations.

    Examples
    ========

    >>> from sympy.abc import x
    >>> from sympy.core.cache import clear_cache
    >>> from sympy.core.parameters import global_parameters as gp

    >>> gp.evaluate
    True
    >>> x+x
    2*x

    >>> log = []
    >>> def f():
    ...     clear_cache()
    ...     gp.evaluate = False
    ...     log.append(x+x)
    ...     clear_cache()
    >>> import threading
    >>> thread = threading.Thread(target=f)
    >>> thread.start()
    >>> thread.join()

    >>> print(log)
    [x + x]

    >>> gp.evaluate
    True
    >>> x+x
    2*x

    References
    ==========

    .. [1] https://docs.python.org/3/library/threading.html

    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __setattr__(self, name, value):
        if getattr(self, name) != value:
            clear_cache()
        return super().__setattr__(name, value)

global_parameters = _global_parameters(evaluate=True, distribute=True)

@contextmanager
def evaluate(x):
    """ Control automatic evaluation

    This context manager controls whether or not all SymPy functions evaluate
    by default.

    Note that much of SymPy expects evaluated expressions.  This functionality
    is experimental and is unlikely to function as intended on large
    expressions.

    Examples
    ========

    >>> from sympy.abc import x
    >>> from sympy.core.parameters import evaluate
    >>> print(x + x)
    2*x
    >>> with evaluate(False):
    ...     print(x + x)
    x + x
    """

    old = global_parameters.evaluate

    try:
        global_parameters.evaluate = x
        yield
    finally:
        global_parameters.evaluate = old


@contextmanager
def distribute(x):
    """ Control automatic distribution of Number over Add

    This context manager controls whether or not Mul distribute Number over
    Add. Plan is to avoid distributing Number over Add in all of sympy. Once
    that is done, this contextmanager will be removed.

    Examples
    ========

    >>> from sympy.abc import x
    >>> from sympy.core.parameters import distribute
    >>> print(2*(x + 1))
    2*x + 2
    >>> with distribute(False):
    ...     print(2*(x + 1))
    2*(x + 1)
    """

    old = global_parameters.distribute

    try:
        global_parameters.distribute = x
        yield
    finally:
        global_parameters.distribute = old
