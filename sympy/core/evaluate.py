from contextlib import contextmanager


global_evaluate = [True]


@contextmanager
def evaluate(x):
    """ Control automatic evaluation

    This context managers controls whether or not all SymPy functions evaluate
    by default.

    Note that much of SymPy expects evaluated expressions.  This functionality
    is experimental and is unlikely to function as intended on large
    expressions.

    Examples
    ========

    >>> from sympy.abc import x
    >>> from sympy.core.evaluate import evaluate
    >>> print(x + x)
    2*x
    >>> with evaluate(False):
    ...     print(x + x)
    x + x
    """

    old = global_evaluate[0]

    global_evaluate[0] = x
    yield
    global_evaluate[0] = old
