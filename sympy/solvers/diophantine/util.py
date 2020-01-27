from sympy.core.compatibility import as_int
from sympy.core.mul import Mul
from sympy.core.numbers import igcd
from sympy.core.power import integer_nthroot
from sympy.functions.elementary.complexes import sign
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.ntheory.factor_ import (
    factorint, multiplicity, perfect_power, divisors)
from sympy.ntheory.generate import nextprime
from sympy.ntheory.primetest import is_square, isprime
from sympy.utilities.misc import filldedent
from sympy.ntheory.residue_ntheory import sqrt_mod


def is_int(i):
    try:
        as_int(i)
        return True
    except ValueError:
        pass


def sorted_tuple(*i):
    return tuple(sorted(i))


def remove_gcd(*x):
    try:
        g = igcd(*x)
    except ValueError:
        fx = list(filter(None, x))
        if len(fx) < 2:
            return x
        g = igcd(*[i.as_content_primitive()[0] for i in fx])
    except TypeError:
        raise TypeError('_remove_gcd(a,b,c) or _remove_gcd(*container)')
    if g == 1:
        return x
    return tuple([i // g for i in x])


def rational_pq(a, b):
    # return `(numer, denom)` for a/b; sign in numer and gcd removed
    return remove_gcd(sign(b) * a, abs(b))


def nint_or_floor(p, q):
    # return nearest int to p/q; in case of tie return floor(p/q)
    w, r = divmod(p, q)
    if abs(r) <= abs(q) // 2:
        return w
    return w + 1


def odd(i):
    return i % 2 != 0


def even(i):
    return i % 2 == 0


def divisible(a, b):
    """
    Returns `True` if ``a`` is divisible by ``b`` and `False` otherwise.
    """
    return not a % b


def length(P, Q, D):
    r"""
    Returns the (length of aperiodic part + length of periodic part) of
    continued fraction representation of `\\frac{P + \sqrt{D}}{Q}`.

    It is important to remember that this does NOT return the length of the
    periodic part but the sum of the lengths of the two parts as mentioned
    above.

    Usage
    =====

    ``length(P, Q, D)``: ``P``, ``Q`` and ``D`` are integers corresponding to
    the continued fraction `\\frac{P + \sqrt{D}}{Q}`.

    Details
    =======

    ``P``, ``D`` and ``Q`` corresponds to P, D and Q in the continued fraction,
    `\\frac{P + \sqrt{D}}{Q}`.

    Examples
    ========

    >>> from sympy.solvers.diophantine.util import length
    >>> length(-2 , 4, 5) # (-2 + sqrt(5))/4
    3
    >>> length(-5, 4, 17) # (-5 + sqrt(17))/4
    4

    See Also
    ========
    sympy.ntheory.continued_fraction.continued_fraction_periodic
    """
    from sympy.ntheory.continued_fraction import continued_fraction_periodic
    v = continued_fraction_periodic(P, Q, D)
    if type(v[-1]) is list:
        rpt = len(v[-1])
        nonrpt = len(v) - 1
    else:
        rpt = 0
        nonrpt = len(v)
    return rpt + nonrpt


def square_factor(a):
    r"""
    Returns an integer `c` s.t. `a = c^2k, \ c,k \in Z`. Here `k` is square
    free. `a` can be given as an integer or a dictionary of factors.

    Examples
    ========

    >>> from sympy.solvers.diophantine.util import square_factor
    >>> square_factor(24)
    2
    >>> square_factor(-36*3)
    6
    >>> square_factor(1)
    1
    >>> square_factor({3: 2, 2: 1, -1: 1})  # -18
    3

    See Also
    ========
    sympy.ntheory.factor_.core
    """
    f = a if isinstance(a, dict) else factorint(a)
    return Mul(*[p ** (e // 2) for p, e in f.items()])


def gaussian_reduce(w, a, b):
    r"""
    Returns a reduced solution `(x, z)` to the congruence
    `X^2 - aZ^2 \equiv 0 \ (mod \ b)` so that `x^2 + |a|z^2` is minimal.

    Details
    =======

    Here ``w`` is a solution of the congruence `x^2 \equiv a \ (mod \ b)`

    References
    ==========

    .. [1] Gaussian lattice Reduction [online]. Available:
           http://home.ie.cuhk.edu.hk/~wkshum/wordpress/?p=404
    .. [2] Efficient Solution of Rational Conices, J. E. Cremona and D. Rusin,
           Mathematics of Computation, Volume 00, Number 0.
    """
    u = (0, 1)
    v = (1, 0)

    if dot(u, v, w, a, b) < 0:
        v = (-v[0], -v[1])

    if norm(u, w, a, b) < norm(v, w, a, b):
        u, v = v, u

    while norm(u, w, a, b) > norm(v, w, a, b):
        k = dot(u, v, w, a, b) // dot(v, v, w, a, b)
        u, v = v, (u[0] - k * v[0], u[1] - k * v[1])

    u, v = v, u

    if dot(u, v, w, a, b) < dot(v, v, w, a, b) / 2 or norm((u[0] - v[0], u[1] - v[1]), w, a, b) > norm(v, w, a, b):
        c = v
    else:
        c = (u[0] - v[0], u[1] - v[1])

    return c[0] * w + b * c[1], c[0]


def dot(u, v, w, a, b):
    r"""
    Returns a special dot product of the vectors `u = (u_{1}, u_{2})` and
    `v = (v_{1}, v_{2})` which is defined in order to reduce solution of
    the congruence equation `X^2 - aZ^2 \equiv 0 \ (mod \ b)`.
    """
    u_1, u_2 = u
    v_1, v_2 = v
    return (w * u_1 + b * u_2) * (w * v_1 + b * v_2) + abs(a) * u_1 * v_1


def norm(u, w, a, b):
    r"""
    Returns the norm of the vector `u = (u_{1}, u_{2})` under the dot product
    defined by `u \cdot v = (wu_{1} + bu_{2})(w*v_{1} + bv_{2}) + |a|*u_{1}*v_{1}`
    where `u = (u_{1}, u_{2})` and `v = (v_{1}, v_{2})`.
    """
    u_1, u_2 = u
    return sqrt(dot((u_1, u_2), (u_1, u_2), w, a, b))



def descent(A, B):
    """
    Returns a non-trivial solution, (x, y, z), to `x^2 = Ay^2 + Bz^2`
    using Lagrange's descent method with lattice-reduction. `A` and `B`
    are assumed to be valid for such a solution to exist.

    This is faster than the normal Lagrange's descent algorithm because
    the Gaussian reduction is used.

    Examples
    ========

    >>> from sympy.solvers.diophantine.util import descent
    >>> descent(3, 1) # x**2 = 3*y**2 + z**2
    (1, 0, 1)

    `(x, y, z) = (1, 0, 1)` is a solution to the above equation.

    >>> descent(41, -113)
    (-16, -3, 1)

    References
    ==========

    .. [1] Efficient Solution of Rational Conices, J. E. Cremona and D. Rusin,
           Mathematics of Computation, Volume 00, Number 0.
    """
    if abs(A) > abs(B):
        x, y, z = descent(B, A)
        return x, z, y

    if B == 1:
        return (1, 0, 1)
    if A == 1:
        return (1, 1, 0)
    if B == -A:
        return (0, 1, 1)
    if B == A:
        x, z, y = descent(-1, A)
        return (A*y, z, x)

    w = sqrt_mod(A, B)
    x_0, z_0 = gaussian_reduce(w, A, B)

    t = (x_0**2 - A*z_0**2) // B
    t_2 = square_factor(t)
    t_1 = t // t_2**2

    x_1, z_1, y_1 = descent(A, t_1)

    return remove_gcd(x_0*x_1 + A*z_0*z_1, z_0*x_1 + x_0*z_1, t_1*t_2*y_1)



def ldescent(A, B):
    """
    Return a non-trivial solution to `w^2 = Ax^2 + By^2` using
    Lagrange's method; return None if there is no such solution.
    .

    Here, `A \\neq 0` and `B \\neq 0` and `A` and `B` are square free. Output a
    tuple `(w_0, x_0, y_0)` which is a solution to the above equation.

    Examples
    ========

    >>> from sympy.solvers.diophantine.util import ldescent
    >>> ldescent(1, 1) # w^2 = x^2 + y^2
    (1, 1, 0)
    >>> ldescent(4, -7) # w^2 = 4x^2 - 7y^2
    (2, -1, 0)

    This means that `x = -1, y = 0` and `w = 2` is a solution to the equation
    `w^2 = 4x^2 - 7y^2`

    >>> ldescent(5, -1) # w^2 = 5x^2 - y^2
    (2, 1, -1)

    References
    ==========

    .. [1] The algorithmic resolution of Diophantine equations, Nigel P. Smart,
           London Mathematical Society Student Texts 41, Cambridge University
           Press, Cambridge, 1998.
    .. [2] Efficient Solution of Rational Conices, J. E. Cremona and D. Rusin,
           [online], Available:
           http://eprints.nottingham.ac.uk/60/1/kvxefz87.pdf
    """
    if abs(A) > abs(B):
        w, y, x = ldescent(B, A)
        return w, x, y

    if A == 1:
        return (1, 1, 0)

    if B == 1:
        return (1, 0, 1)

    if B == -1:  # and A == -1
        return

    r = sqrt_mod(A, B)

    Q = (r**2 - A) // B

    if Q == 0:
        B_0 = 1
        d = 0
    else:
        div = divisors(Q)
        B_0 = None

        for i in div:
            sQ, _exact = integer_nthroot(abs(Q) // i, 2)
            if _exact:
                B_0, d = sign(Q)*i, sQ
                break

    if B_0 is not None:
        W, X, Y = ldescent(A, B_0)
        return remove_gcd((-A*X + r*W), (r*X - W), Y*(B_0*d))


## Functions below this comment can be more suitably grouped under
## an Additive number theory module rather than the Diophantine
## equation module.


def partition(n, k=None, zeros=False):
    """
    Returns a generator that can be used to generate partitions of an integer
    `n`.

    A partition of `n` is a set of positive integers which add up to `n`. For
    example, partitions of 3 are 3, 1 + 2, 1 + 1 + 1. A partition is returned
    as a tuple. If ``k`` equals None, then all possible partitions are returned
    irrespective of their size, otherwise only the partitions of size ``k`` are
    returned. If the ``zero`` parameter is set to True then a suitable
    number of zeros are added at the end of every partition of size less than
    ``k``.

    ``zero`` parameter is considered only if ``k`` is not None. When the
    partitions are over, the last `next()` call throws the ``StopIteration``
    exception, so this function should always be used inside a try - except
    block.

    Details
    =======

    ``partition(n, k)``: Here ``n`` is a positive integer and ``k`` is the size
    of the partition which is also positive integer.

    Examples
    ========

    >>> from sympy.solvers.diophantine.util import partition
    >>> f = partition(5)
    >>> next(f)
    (1, 1, 1, 1, 1)
    >>> next(f)
    (1, 1, 1, 2)
    >>> g = partition(5, 3)
    >>> next(g)
    (1, 1, 3)
    >>> next(g)
    (1, 2, 2)
    >>> g = partition(5, 3, zeros=True)
    >>> next(g)
    (0, 0, 5)

    """
    from sympy.utilities.iterables import ordered_partitions
    if not zeros or k is None:
        for i in ordered_partitions(n, k):
            yield tuple(i)
    else:
        for m in range(1, k + 1):
            for i in ordered_partitions(n, m):
                i = tuple(i)
                yield (0,) * (k - len(i)) + i


def prime_as_sum_of_two_squares(p):
    """
    Represent a prime `p` as a unique sum of two squares; this can
    only be done if the prime is congruent to 1 mod 4.

    Examples
    ========

    >>> from sympy.solvers.diophantine.util import prime_as_sum_of_two_squares
    >>> prime_as_sum_of_two_squares(7)  # can't be done
    >>> prime_as_sum_of_two_squares(5)
    (1, 2)

    Reference
    =========

    .. [1] Representing a number as a sum of four squares, [online],
        Available: http://schorn.ch/lagrange.html

    See Also
    ========
    sum_of_squares()
    """
    if not p % 4 == 1:
        return

    if p % 8 == 5:
        b = 2
    else:
        b = 3

        while pow(b, (p - 1) // 2, p) == 1:
            b = nextprime(b)

    b = pow(b, (p - 1) // 4, p)
    a = p

    while b ** 2 > p:
        a, b = b, a % b

    return (int(a % b), int(b))  # convert from long


def sum_of_three_squares(n):
    r"""
    Returns a 3-tuple `(a, b, c)` such that `a^2 + b^2 + c^2 = n` and
    `a, b, c \geq 0`.

    Returns None if `n = 4^a(8m + 7)` for some `a, m \in Z`. See
    [1]_ for more details.

    Usage
    =====

    ``sum_of_three_squares(n)``: Here ``n`` is a non-negative integer.

    Examples
    ========

    >>> from sympy.solvers.diophantine.util import sum_of_three_squares
    >>> sum_of_three_squares(44542)
    (18, 37, 207)

    References
    ==========

    .. [1] Representing a number as a sum of three squares, [online],
        Available: http://schorn.ch/lagrange.html

    See Also
    ========
    sum_of_squares()
    """
    special = {1: (1, 0, 0), 2: (1, 1, 0), 3: (1, 1, 1), 10: (1, 3, 0), 34: (3, 3, 4), 58: (3, 7, 0),
               85: (6, 7, 0), 130: (3, 11, 0), 214: (3, 6, 13), 226: (8, 9, 9), 370: (8, 9, 15),
               526: (6, 7, 21), 706: (15, 15, 16), 730: (1, 27, 0), 1414: (6, 17, 33), 1906: (13, 21, 36),
               2986: (21, 32, 39), 9634: (56, 57, 57)}

    v = 0

    if n == 0:
        return (0, 0, 0)

    v = multiplicity(4, n)
    n //= 4 ** v

    if n % 8 == 7:
        return

    if n in special.keys():
        x, y, z = special[n]
        return sorted_tuple(2 ** v * x, 2 ** v * y, 2 ** v * z)

    s, _exact = integer_nthroot(n, 2)

    if _exact:
        return (2 ** v * s, 0, 0)

    x = None

    if n % 8 == 3:
        s = s if odd(s) else s - 1

        for x in range(s, -1, -2):
            N = (n - x ** 2) // 2
            if isprime(N):
                y, z = prime_as_sum_of_two_squares(N)
                return sorted_tuple(2 ** v * x, 2 ** v * (y + z), 2 ** v * abs(y - z))
        return

    if n % 8 == 2 or n % 8 == 6:
        s = s if odd(s) else s - 1
    else:
        s = s - 1 if odd(s) else s

    for x in range(s, -1, -2):
        N = n - x ** 2
        if isprime(N):
            y, z = prime_as_sum_of_two_squares(N)
            return sorted_tuple(2 ** v * x, 2 ** v * y, 2 ** v * z)


def sum_of_four_squares(n):
    r"""
    Returns a 4-tuple `(a, b, c, d)` such that `a^2 + b^2 + c^2 + d^2 = n`.

    Here `a, b, c, d \geq 0`.

    Usage
    =====

    ``sum_of_four_squares(n)``: Here ``n`` is a non-negative integer.

    Examples
    ========

    >>> from sympy.solvers.diophantine.util import sum_of_four_squares
    >>> sum_of_four_squares(3456)
    (8, 8, 32, 48)
    >>> sum_of_four_squares(1294585930293)
    (0, 1234, 2161, 1137796)

    References
    ==========

    .. [1] Representing a number as a sum of four squares, [online],
        Available: http://schorn.ch/lagrange.html

    See Also
    ========
    sum_of_squares()
    """
    if n == 0:
        return (0, 0, 0, 0)

    v = multiplicity(4, n)
    n //= 4 ** v

    if n % 8 == 7:
        d = 2
        n = n - 4
    elif n % 8 == 6 or n % 8 == 2:
        d = 1
        n = n - 1
    else:
        d = 0

    x, y, z = sum_of_three_squares(n)

    return sorted_tuple(2 ** v * d, 2 ** v * x, 2 ** v * y, 2 ** v * z)


def power_representation(n, p, k, zeros=False):
    r"""
    Returns a generator for finding k-tuples of integers,
    `(n_{1}, n_{2}, . . . n_{k})`, such that
    `n = n_{1}^p + n_{2}^p + . . . n_{k}^p`.

    Usage
    =====

    ``power_representation(n, p, k, zeros)``: Represent non-negative number
    ``n`` as a sum of ``k`` ``p``\ th powers. If ``zeros`` is true, then the
    solutions is allowed to contain zeros.

    Examples
    ========

    >>> from sympy.solvers.diophantine.util import power_representation

    Represent 1729 as a sum of two cubes:

    >>> f = power_representation(1729, 3, 2)
    >>> next(f)
    (9, 10)
    >>> next(f)
    (1, 12)

    If the flag `zeros` is True, the solution may contain tuples with
    zeros; any such solutions will be generated after the solutions
    without zeros:

    >>> list(power_representation(125, 2, 3, zeros=True))
    [(5, 6, 8), (3, 4, 10), (0, 5, 10), (0, 2, 11)]

    For even `p` the `permute_sign` function can be used to get all
    signed values:

    >>> from sympy.utilities.iterables import permute_signs
    >>> list(permute_signs((1, 12)))
    [(1, 12), (-1, 12), (1, -12), (-1, -12)]

    All possible signed permutations can also be obtained:

    >>> from sympy.utilities.iterables import signed_permutations
    >>> list(signed_permutations((1, 12)))
    [(1, 12), (-1, 12), (1, -12), (-1, -12), (12, 1), (-12, 1), (12, -1), (-12, -1)]
    """
    n, p, k = [as_int(i) for i in (n, p, k)]

    if n < 0:
        if p % 2:
            for t in power_representation(-n, p, k, zeros):
                yield tuple(-i for i in t)
        return

    if p < 1 or k < 1:
        raise ValueError(filldedent('''
    Expecting positive integers for `(p, k)`, but got `(%s, %s)`'''
                                    % (p, k)))

    if n == 0:
        if zeros:
            yield (0,) * k
        return

    if k == 1:
        if p == 1:
            yield (n,)
        else:
            be = perfect_power(n)
            if be:
                b, e = be
                d, r = divmod(e, p)
                if not r:
                    yield (b ** d,)
        return

    if p == 1:
        for t in partition(n, k, zeros=zeros):
            yield t
        return

    if p == 2:
        feasible = can_do_sum_of_squares(n, k)
        if not feasible:
            return
        if not zeros and n > 33 and k >= 5 and k <= n and n - k in (
            13, 10, 7, 5, 4, 2, 1):
            '''Todd G. Will, "When Is n^2 a Sum of k Squares?", [online].
                Available: https://www.maa.org/sites/default/files/Will-MMz-201037918.pdf'''
            return
        if feasible is not True:  # it's prime and k == 2
            yield prime_as_sum_of_two_squares(n)
            return

    if k == 2 and p > 2:
        be = perfect_power(n)
        if be and be[1] % p == 0:
            return  # Fermat: a**n + b**n = c**n has no solution for n > 2

    if n >= k:
        a = integer_nthroot(n - (k - 1), p)[0]
        for t in pow_rep_recursive(a, k, n, [], p):
            yield tuple(reversed(t))

    if zeros:
        a = integer_nthroot(n, p)[0]
        for i in range(1, k):
            for t in pow_rep_recursive(a, i, n, [], p):
                yield tuple(reversed(t + (0,) * (k - i)))


def pow_rep_recursive(n_i, k, n_remaining, terms, p):
    if k == 0 and n_remaining == 0:
        yield tuple(terms)
    else:
        if n_i >= 1 and k > 0:
            for t in pow_rep_recursive(n_i - 1, k, n_remaining, terms, p):
                yield t
            residual = n_remaining - pow(n_i, p)
            if residual >= 0:
                for t in pow_rep_recursive(n_i, k - 1, residual, terms + [n_i], p):
                    yield t


def sum_of_squares(n, k, zeros=False):
    """Return a generator that yields the k-tuples of nonnegative
    values, the squares of which sum to n. If zeros is False (default)
    then the solution will not contain zeros. The nonnegative
    elements of a tuple are sorted.

    * If k == 1 and n is square, (n,) is returned.

    * If k == 2 then n can only be written as a sum of squares if
      every prime in the factorization of n that has the form
      4*k + 3 has an even multiplicity. If n is prime then
      it can only be written as a sum of two squares if it is
      in the form 4*k + 1.

    * if k == 3 then n can be written as a sum of squares if it does
      not have the form 4**m*(8*k + 7).

    * all integers can be written as the sum of 4 squares.

    * if k > 4 then n can be partitioned and each partition can
      be written as a sum of 4 squares; if n is not evenly divisible
      by 4 then n can be written as a sum of squares only if the
      an additional partition can be written as sum of squares.
      For example, if k = 6 then n is partitioned into two parts,
      the first being written as a sum of 4 squares and the second
      being written as a sum of 2 squares -- which can only be
      done if the condition above for k = 2 can be met, so this will
      automatically reject certain partitions of n.

    Examples
    ========

    >>> from sympy.solvers.diophantine.util import sum_of_squares
    >>> list(sum_of_squares(25, 2))
    [(3, 4)]
    >>> list(sum_of_squares(25, 2, True))
    [(3, 4), (0, 5)]
    >>> list(sum_of_squares(25, 4))
    [(1, 2, 2, 4)]

    See Also
    ========
    sympy.utilities.iterables.signed_permutations
    """
    for t in power_representation(n, 2, k, zeros):
        yield t


def can_do_sum_of_squares(n, k):
    """Return True if n can be written as the sum of k squares,
    False if it can't, or 1 if k == 2 and n is prime (in which
    case it *can* be written as a sum of two squares). A False
    is returned only if it can't be written as k-squares, even
    if 0s are allowed.
    """
    if k < 1:
        return False
    if n < 0:
        return False
    if n == 0:
        return True
    if k == 1:
        return is_square(n)
    if k == 2:
        if n in (1, 2):
            return True
        if isprime(n):
            if n % 4 == 1:
                return 1  # signal that it was prime
            return False
        else:
            f = factorint(n)
            for p, m in f.items():
                # we can proceed iff no prime factor in the form 4*k + 3
                # has an odd multiplicity
                if (p % 4 == 3) and m % 2:
                    return False
            return True
    if k == 3:
        if (n // 4 ** multiplicity(4, n)) % 8 == 7:
            return False
    # every number can be written as a sum of 4 squares; for k > 4 partitions
    # can be 0
    return True
