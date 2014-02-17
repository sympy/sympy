from contextlib import contextmanager


global_evaluate = [True]


@contextmanager
def evaluate(x):
    """ Switch automatic evaluation on or off

    Examples
    ========

    >>> from sympy.abc import x
    >>> from sympy.core.operations import evaluate
    >>> print x + x
    2*x
    >>> with evaluate(False):
    ...     print x + x
    x + x
    """

    old = global_evaluate[0]

    global_evaluate[0] = x
    yield
    global_evaluate[0] = old
