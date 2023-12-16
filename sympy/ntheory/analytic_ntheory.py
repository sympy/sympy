from .factor_ import factorint
from sympy.utilities.memoization import recurrence_memo
from sympy.utilities.misc import as_int


@recurrence_memo([0, 1])
def _ramanujan_tau(n, prev):
    r""" Calculate Ramanujan tau function `\tau(n)`.

    Explanation
    ===========

    This is a low-level helper for ``ramanujan_tau``, for internal use.
    The following equation is used to obtain `\tau(n)`.

    .. math ::
        (n-1)\tau(n) = \sum_{m=1}^{b_n} (-1)^{m+1} (n-1-\frac{9}{2}m(m+1))\tau(n-\frac{1}{2}m(m+1))

    Therefore, the results of calculations below `n` are memoized and used.

    Parameters
    ==========

    n : int
        positive integer

    Returns
    =======

    int : `\tau(n)`

    Examples
    ========

    >>> from sympy.ntheory.analytic_ntheory import _ramanujan_tau
    >>> _ramanujan_tau(2)
    -24

    References
    ==========

    .. [1] https://mathworld.wolfram.com/TauFunction.html

    """
    t = 0
    u = 0
    nm1 = n - 1
    # The end of the loop is dummy.
    # In fact, exit the for statement with the condition `n <= u`.
    for m in range(1, n):
        u += m
        if n <= u:
            break
        t += (1 if m % 2 else -1)*(2*m + 1)*(nm1 - 9*u)*prev[-u]
    return t // nm1


def ramanujan_tau(n):
    r""" Calculate Ramanujan tau function `\tau(n)`.

    Explanation
    ===========

    `\tau(n)` is a multiplicative function.
    That is, if `\gcd(a,b)=1` then `\tau(ab)=\tau(a)\tau(b)`.
    Furthermore, when `p` is a prime number and `r>1`, the following recurrence relation holds.

    .. math ::
        \tau(p^r) = \tau(p)\tau(p^{r-1}) - p^{11}\tau(p^{r-2})

    Therefore, prime factorize `n` to obtain `\tau(p)` and then `\tau(p^r)`.
    The product of these is `\tau(n)`.

    Parameters
    ==========

    n : int
        positive integer

    Returns
    =======

    int : `\tau(n)`

    Raises
    ======

    ValueError
        If n is not a positive integer

    Examples
    ========

    >>> from sympy.ntheory.analytic_ntheory import ramanujan_tau
    >>> ramanujan_tau(10)
    -115920

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Ramanujan_tau_function
    .. [2] https://oeis.org/A000594

    """
    n = as_int(n)
    if n < 1:
        raise ValueError("n should be a positive integer")
    if n < _ramanujan_tau.cache_length():
        return _ramanujan_tau(n)
    result = 1
    for p, e in factorint(n).items():
        t1 = tau_p = _ramanujan_tau(p)
        t2 = 1 # tau(p^0)
        for _ in range(e - 1):
            t1, t2 = tau_p*t1 - p**11*t2, t1
        result *= t1
    return result
