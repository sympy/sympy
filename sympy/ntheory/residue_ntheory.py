from sympy.core.numbers import igcd
from primetest import isprime
from factor_ import factorint


def totient_(n):
    """returns the number of integers less than n
    and relatively prime to n"""
    if n < 1:
        raise ValueError("n must be a positive integer")
    tot = 0
    for x in xrange(1, n):
        if igcd(x, n) == 1:
            tot += 1
    return tot


def n_order(a, n):
    """ returns the order of a modulo n
    Order of a modulo n is the smallest integer
    k such that a^k leaves a remainder of 1 with n.
    """
    if igcd(a, n) != 1:
        raise ValueError("The two numbers should be relatively prime")
    group_order = totient_(n)
    factors = factorint(group_order)
    order = 1
    if a > n:
        a = a % n
    for p, e in factors.iteritems():
        exponent = group_order
        for f in xrange(0, e + 1):
            if (a ** (exponent)) % n != 1:
                order *= p ** (e - f + 1)
                break
            exponent = exponent // p
    return order


def is_primitive_root(a, p):
    """
    returns True if a is a primitive root of p
    """
    if igcd(a, p) != 1:
        raise ValueError("The two numbers should be relatively prime")
    if a > p:
        a = a % p
    if n_order(a, p) == totient_(p):
        return True
    else:
        return False


def is_quad_residue(a, p):
    """
    returns True if a is a quadratic residue of p
    p should be a prime and a should be relatively
    prime to p
    """
    if not isprime(p) or p == 2:
        raise ValueError("p should be an odd prime")
    if igcd(a, p) != 1:
        raise ValueError("The two numbers should be relatively prime")
    if a > p:
        a = a % p

    def square_and_multiply(a, n, p):
        if n == 0:
            return 1
        elif n == 1:
            return a
        elif n % 2 == 1:
            return ((square_and_multiply(a, n // 2, p) ** 2) * a) % p
        else:
            return (square_and_multiply(a, n // 2, p) ** 2) % p

    return (square_and_multiply(a, (p - 1) // 2, p) % p) == 1


def legendre_symbol(a, p):
    """
    return 1 if a is a quadratic residue of p
    else return -1
    p should be an odd prime by definition
    """
    if not isprime(p) or p == 2:
        raise ValueError("p should be an odd prime")
    if igcd(a, p) != 1:
        raise ValueError("The two numbers should be relatively prime")
    if a > p:
        a = a % p
    if is_quad_residue(a, p):
        return 1
    else:
        return -1
