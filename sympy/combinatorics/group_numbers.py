from sympy.ntheory.factor_ import factorint
from sympy.utilities.misc import as_int


def _is_nilpotent_number(factors: dict) -> bool:
    """ Check whether `n` is a nilpotent number.
    Note that ``factors`` is a prime factorization of `n`.

    This is a low-level helper for ``is_nilpotent_number``, for internal use.
    """
    for p in factors.keys():
        for q, e in factors.items():
            # We want to calculate
            # any(pow(q, k, p) == 1 for k in range(1, e + 1))
            m = 1
            for _ in range(e):
                m = m*q % p
                if m == 1:
                    return False
    return True


def is_nilpotent_number(n) -> bool:
    """
    Check whether `n` is a nilpotent number. A number `n` is said to be
    nilpotent if and only if every finite group of order `n` is nilpotent.
    For more information see [1]_.

    Examples
    ========

    >>> from sympy.combinatorics.group_numbers import is_nilpotent_number
    >>> from sympy import randprime
    >>> is_nilpotent_number(21)
    False
    >>> is_nilpotent_number(randprime(1, 30)**12)
    True

    References
    ==========

    .. [1] Pakianathan, J., Shankar, K., Nilpotent Numbers,
           The American Mathematical Monthly, 107(7), 631-634.
    .. [2] https://oeis.org/A056867

    """
    n = as_int(n)
    if n <= 0:
        raise ValueError("n must be a positive integer, not %i" % n)
    return _is_nilpotent_number(factorint(n))


def is_abelian_number(n) -> bool:
    """
    Check whether `n` is an abelian number. A number `n` is said to be abelian
    if and only if every finite group of order `n` is abelian. For more
    information see [1]_.

    Examples
    ========

    >>> from sympy.combinatorics.group_numbers import is_abelian_number
    >>> from sympy import randprime
    >>> is_abelian_number(4)
    True
    >>> is_abelian_number(randprime(1, 2000)**2)
    True
    >>> is_abelian_number(60)
    False

    References
    ==========

    .. [1] Pakianathan, J., Shankar, K., Nilpotent Numbers,
           The American Mathematical Monthly, 107(7), 631-634.
    .. [2] https://oeis.org/A051532

    """
    n = as_int(n)
    if n <= 0:
        raise ValueError("n must be a positive integer, not %i" % n)
    factors = factorint(n)
    return all(e < 3 for e in factors.values()) and _is_nilpotent_number(factors)


def is_cyclic_number(n) -> bool:
    """
    Check whether `n` is a cyclic number. A number `n` is said to be cyclic
    if and only if every finite group of order `n` is cyclic. For more
    information see [1]_.

    Examples
    ========

    >>> from sympy.combinatorics.group_numbers import is_cyclic_number
    >>> from sympy import randprime
    >>> is_cyclic_number(15)
    True
    >>> is_cyclic_number(randprime(1, 2000)**2)
    False
    >>> is_cyclic_number(4)
    False

    References
    ==========

    .. [1] Pakianathan, J., Shankar, K., Nilpotent Numbers,
           The American Mathematical Monthly, 107(7), 631-634.
    .. [2] https://oeis.org/A003277

    """
    n = as_int(n)
    if n <= 0:
        raise ValueError("n must be a positive integer, not %i" % n)
    factors = factorint(n)
    return all(e == 1 for e in factors.values()) and _is_nilpotent_number(factors)
