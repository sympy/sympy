from sympy.core import Integer, Pow, Mod
from sympy import factorint


def is_nilpotent_number(n):
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

    .. [1] Pakianathan, J., Shankar, K., *Nilpotent Numbers*,
            The American Mathematical Monthly, 107(7), 631-634.


    """
    if n <= 0 or int(n) != n:
        raise ValueError("n must be a positive integer, not %i" % n)

    n = Integer(n)
    prime_factors = list(factorint(n).items())
    is_nilpotent = True
    for p_j, a_j in prime_factors:
        for p_i, a_i in prime_factors:
            if any([Mod(Pow(p_i, k), p_j) == 1 for k in range(1, a_i + 1)]):
                is_nilpotent = False
                break
        if not is_nilpotent:
            break

    return is_nilpotent


def is_abelian_number(n):
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

    .. [1] Pakianathan, J., Shankar, K., *Nilpotent Numbers*,
            The American Mathematical Monthly, 107(7), 631-634.


    """
    if n <= 0 or int(n) != n:
        raise ValueError("n must be a positive integer, not %i" % n)

    n = Integer(n)
    if not is_nilpotent_number(n):
        return False

    prime_factors = list(factorint(n).items())
    is_abelian = all(a_i < 3 for p_i, a_i in prime_factors)
    return is_abelian


def is_cyclic_number(n):
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

    .. [1] Pakianathan, J., Shankar, K., *Nilpotent Numbers*,
            The American Mathematical Monthly, 107(7), 631-634.

    """
    if n <= 0 or int(n) != n:
        raise ValueError("n must be a positive integer, not %i" % n)

    n = Integer(n)
    if not is_nilpotent_number(n):
        return False

    prime_factors = list(factorint(n).items())
    is_cyclic = all(a_i < 2 for p_i, a_i in prime_factors)
    return is_cyclic
