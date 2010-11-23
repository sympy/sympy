from sympy.core.numbers import igcd
from primetest import isprime

def totient_(n):
    """returns the number of integers less than n
    and relatively prime to n"""
    if n < 1:
        raise ValueError("n must be a positive integer")
    tot=0
    for x in xrange(1,n):
        if igcd(x,n)==1:
            tot+=1
    return tot

def n_order(a,n):
    """ returns the order of a modulo n
    Order of a modulo n is the smallest integer
    k such that a^k leaves a remainder of 1 with n.
    """
    assert igcd(a,n)==1
    if a>n : a=a%n
    for x in xrange(1,totient_(n)+1):
        if (a**x)%n==1:
            return x

def is_primitive_root(a,p):
    """
    returns True if a is a primitive root of p
    """
    assert igcd(a,p) == 1,"The two numbers should be relatively prime"
    if a>p:
        a=a%p
    if n_order(a,p)==totient_(p):
        return True
    else:
        return False

def is_quad_residue(a,p):
    """
    returns True if a is a quadratic residue of p
    p should be a prime and a should be relatively
    prime to p
    """
    assert isprime(p) and p!=2,"p should be an odd prime"
    assert igcd(a,p)==1,"The two numbers should be relatively prime"
    if a>p:
        a=a%p
    rem=(a**((p-1)//2))%p    # a^(p-1 / 2) % p
    if rem==1: return True
    else : return False

def legendre_symbol(a,p):
    """
    return 1 if a is a quadratic residue of p
    else return -1
    p should be an odd prime by definition
    """
    assert isprime(p) and p!=2,"p should be an odd prime"
    assert igcd(a,p)==1,"The two numbers should be relatively prime"
    if a>p:
        a=a%p
    if is_quad_residue(a,p)==True: return 1
    else : return -1
