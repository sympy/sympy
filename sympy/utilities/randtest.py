""" Helpers for randomized testing """

from __future__ import print_function, division

from random import uniform
import random

from sympy import I, nsimplify, Tuple, Symbol, sympify
from sympy.core.compatibility import is_sequence, as_int


def random_complex_number(a=2, b=-1, c=3, d=1, rational=False):
    """
    Return a random complex number.

    To reduce chance of hitting branch cuts or anything, we guarantee
    b <= Im z <= d, a <= Re z <= c
    """
    A, B = uniform(a, c), uniform(b, d)
    if not rational:
        return A + I*B
    return nsimplify(A, rational=True) + I*nsimplify(B, rational=True)


def comp(z1, z2, tol):
    """Return a bool indicating whether the error between z1 and z2 is <= tol.

    If z2 is non-zero and ``|z1| > 1`` the error is normalized by ``|z1|``, so
    if you want the absolute error, call this as ``comp(z1 - z2, 0, tol)``.
    """
    if not z1:
        z1, z2 = z2, z1
    if not z1:
        return True
    diff = abs(z1 - z2)
    az1 = abs(z1)
    if z2 and az1 > 1:
        return diff/az1 <= tol
    else:
        return diff <= tol


def numerically_zero(f, number=10, prec=10, a=2, b=-1, c=3, d=1, maxfail=100, maxprec=100):
    """
    Test numerically if ``|f| < 10**-prec`` by replacing all symbols (``number``
    times) with random complex numbers and evaluating the result. If a value
    is ever computed that fails the test, False will be returned; otherwise
    None will be returned. Failures to compute a finite result will not be
    counted and will be allowed ``maxfail`` times before an error message
    is raised.

    If you are not sure what precision to use you can set ``prec=None`` and
    a numerical estimate of the error will be made from two numerical
    estimates.

    Examples
    ========

    >>> from sympy import sin, cos, pi, S, Add
    >>> from sympy.abc import x
    >>> from sympy.utilities.randtest import numerically_zero as is_zero
    >>> is_zero(sin(x)**2 + cos(x)**2)
    False
    >>> is_zero(sin(x)**2 + cos(x)**2 - 1)

    Better detection of nonzero values is obtained by using larger values
    of ``prec``:

    >>> eps = lambda e: Add(pi, 10**S(e), -pi, evaluate=False)
    >>> is_zero(eps(-11))
    >>> is_zero(eps(-11), prec=12)
    False

    See Also
    ========
    random_complex_number
    """
    from sympy.core.function import _coeff_isneg
    tol = sympify(10)**-prec if prec else None
    f = sympify(f)
    if not f.args:
        if f.is_number and abs(f) <= (tol or 0):
            return True
        return False
    s = f.free_symbols
    did = 0
    for i in range(number + maxfail):
        reps = dict(zip(s,
            [random_complex_number(a, b, c, d, rational=True) for i in s]))
        if tol:
            check = f.xreplace(reps).n(2, maxn=prec)
            if not check.is_finite == True:
                # reject attempts where nan or infinite quantites were calculated
                continue
            if abs(check) > tol:
                return False
        else:
            M = 10
            n = f.xreplace(reps)
            n2 = n.n(2, maxn=M)
            r, i = n2.as_real_imag()
            # if neither real not imaginary is computed with precision, then
            # retry at maximum working precision
            if r._prec == 1 and i._prec == 1 and maxprec > 10:
                M = maxprec
                n2 = n.n(2, maxn=M)
            n4 = n.n(4, maxn=M)
            err = abs(n4 - n2)
            if abs(n4) > err:
                return False
        did += 1
        if not s or did == number:
            break
    else:
        raise NotImplementedError('unable to find good values for calculation')


def verify_numerically(f, g, z=None, tol=1.0e-6, a=2, b=-1, c=3, d=1):
    """
    Test numerically that f and g agree when evaluated in the argument z.

    If z is None, all symbols will be tested. This routine does not test
    whether there are Floats present with precision higher than 15 digits
    so if there are, your results may not be what you expect due to round-
    off errors.

    Examples
    ========

    >>> from sympy import sin, cos
    >>> from sympy.abc import x
    >>> from sympy.utilities.randtest import verify_numerically as tn
    >>> tn(sin(x)**2 + cos(x)**2, 1, x)
    True
    """
    f, g, z = Tuple(f, g, z)
    z = [z] if isinstance(z, Symbol) else (f.free_symbols | g.free_symbols)
    reps = list(zip(z, [random_complex_number(a, b, c, d) for zi in z]))
    z1 = f.subs(reps).n()
    z2 = g.subs(reps).n()
    return comp(z1, z2, tol)


def test_derivative_numerically(f, z, tol=1.0e-6, a=2, b=-1, c=3, d=1):
    """
    Test numerically that the symbolically computed derivative of f
    with respect to z is correct.

    This routine does not test whether there are Floats present with
    precision higher than 15 digits so if there are, your results may
    not be what you expect due to round-off errors.

    Examples
    ========

    >>> from sympy import sin
    >>> from sympy.abc import x
    >>> from sympy.utilities.randtest import test_derivative_numerically as td
    >>> td(sin(x), x)
    True
    """
    from sympy.core.function import Derivative
    z0 = random_complex_number(a, b, c, d)
    f1 = f.diff(z).subs(z, z0)
    f2 = Derivative(f, z).doit_numerically(z0)
    return comp(f1.n(), f2.n(), tol)

def _randrange(seed=None):
    """Return a randrange generator. ``seed`` can be
        o None - return randomly seeded generator
        o int - return a generator seeded with the int
        o list - the values to be returned will be taken from the list
          in the order given; the provided list is not modified.

    Examples
    ========

    >>> from sympy.utilities.randtest import _randrange
    >>> rr = _randrange()
    >>> rr(1000) # doctest: +SKIP
    999
    >>> rr = _randrange(3)
    >>> rr(1000) # doctest: +SKIP
    238
    >>> rr = _randrange([0, 5, 1, 3, 4])
    >>> rr(3), rr(3)
    (0, 1)
    """
    if seed is None:
        return random.randrange
    elif isinstance(seed, int):
        return random.Random(seed).randrange
    elif is_sequence(seed):
        seed = list(seed)  # make a copy
        seed.reverse()

        def give(a, b=None, seq=seed):
            if b is None:
                a, b = 0, a
            a, b = as_int(a), as_int(b)
            w = b - a
            if w < 1:
                raise ValueError('_randrange got empty range')
            try:
                x = seq.pop()
            except AttributeError:
                raise ValueError('_randrange expects a list-like sequence')
            except IndexError:
                raise ValueError('_randrange sequence was too short')
            if a <= x < b:
                return x
            else:
                return give(a, b, seq)
        return give
    else:
        raise ValueError('_randrange got an unexpected seed')


def _randint(seed=None):
    """Return a randint generator. ``seed`` can be
        o None - return randomly seeded generator
        o int - return a generator seeded with the int
        o list - the values to be returned will be taken from the list
          in the order given; the provided list is not modified.

    Examples
    ========

    >>> from sympy.utilities.randtest import _randint
    >>> ri = _randint()
    >>> ri(1, 1000) # doctest: +SKIP
    999
    >>> ri = _randint(3)
    >>> ri(1, 1000) # doctest: +SKIP
    238
    >>> ri = _randint([0, 5, 1, 2, 4])
    >>> ri(1, 3), ri(1, 3)
    (1, 2)
    """
    if seed is None:
        return random.randint
    elif isinstance(seed, int):
        return random.Random(seed).randint
    elif is_sequence(seed):
        seed = list(seed)  # make a copy
        seed.reverse()

        def give(a, b, seq=seed):
            a, b = as_int(a), as_int(b)
            w = b - a
            if w < 0:
                raise ValueError('_randint got empty range')
            try:
                x = seq.pop()
            except AttributeError:
                raise ValueError('_randint expects a list-like sequence')
            except IndexError:
                raise ValueError('_randint sequence was too short')
            if a <= x <= b:
                return x
            else:
                return give(a, b, seq)
        return give
    else:
        raise ValueError('_randint got an unexpected seed')
