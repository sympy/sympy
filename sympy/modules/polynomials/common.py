"""Simple helper functions common to several algorithms"""

from sympy.core import *
from sympy.core.basic import S # Use Singleton comparisons.
from sympy.core import numbers # Need numbers.gcd
from sympy.core.numbers import NumberSymbol, ImaginaryUnit # To look for numbers
from sympy.modules.utilities import *

from sympy.modules.polynomials.base import PolynomialException, coeff_rings

def reverse(t):
    """Return a tuple with reversed order"""
    return tuple([t[i] for i in range(len(t)-1, 0, -1)])

def term_cmp(a, b, order):
    if order == 'lex':
        return cmp(a[1:], b[1:])
    elif order == 'grlex':
        return cmp((sum(a[1:]),) + a[1:], (sum(b[1:]),) + b[1:])
    elif order == 'grevlex':
        return cmp((sum(a[1:]),) + reverse(map(lambda l:-l, a[1:])),
                   (sum(b[1:]),) + reverse(map(lambda l:-l, b[1:])))
    elif order == '1-el':
        return cmp((a[1], sum(a[2:])) + reverse(map(lambda l:-l,a[2:])),
                   (b[1], sum(b[2:])) + reverse(map(lambda l:-l,b[2:])))
    else:
        raise PolynomialException(str(order) + 'is not an implemented order.')

def term_mult(a, b):
    """Returns a term that represents the multiplication of a and b.

    a and b are assumed to be terms of coefficient lists of
    Polynomials of same the variables.
    """
    return ((a[0]*b[0]).expand(),) \
           + tuple(map(lambda (x,y): x+y, zip(a[1:], b[1:])))

def term_div(a, b):
    """Returns a term that represents the division of a by b.

    a and b are assumed to be terms of coefficient lists of
    Polynomials of same the variables. Divisibility is not tested.
    """
    return ((a[0]/b[0]).expand(),) \
           + tuple(map(lambda (x,y): x-y, zip(a[1:], b[1:])))

def term_is_mult(a, b):
    """Return True if a is a multiple of b

    a and b are assumed to be terms of coefficient lists of
    Polynomials of same the variables."""
    return all([x.is_nonnegative for x in term_div(a, b)[1:]])

def term_lcm(a, b):
    # TODO: Compute lcm oder product of coefficients?
    return (S.One,) + tuple([max(aa, bb)
                             for aa, bb in zip(a[1:], b[1:])])

def merge_var(*a):
    """Return a sorted list of the symbols in the arguments"""
    result = []
    for var in a:
        for sym in var:
            if not sym in result:
                result.append(sym)
    result.sort()
    return result

def coeff_ring(atom):
    """Determine the coefficient ring of some atom, or some list of atoms.
    """
    if isinstance(atom, (Number, NumberSymbol, ImaginaryUnit)) \
        or (isinstance(atom, (Add, Mul))
            and all(map(lambda a:isinstance(a, (Number, NumberSymbol,
                                                ImaginaryUnit)), atom[:]))) \
        or (isinstance(atom, Pow) and isinstance(atom.base, (Number,
                                                             NumberSymbol,
                                                             ImaginaryUnit))
            and isinstance(atom.exp, (Number, NumberSymbol, ImaginaryUnit))):
        if atom.is_integer:
            return 'int'
        elif atom.is_rational:
            return 'rat'
        elif atom.is_real:
            return 'real'
        else:
            return 'cplx'
    elif isinstance(atom, list):
        # Get the coefficient ring of each atom and look for the worst case.
        result = 'int'
        for a in atom:
            cr = coeff_ring(a)
            assert type(cr) == str
            if coeff_rings.index(cr) > coeff_rings.index(result):
                result = cr
        return result
    else:
        return 'sym'

def get_numbers(atom):
    # TODO: Merge with coeff_ring, without recursion!
    result = []
    if isinstance(atom, (Number, NumberSymbol, ImaginaryUnit)) \
        or (isinstance(atom, (Add, Mul))
            and all(map(lambda a:isinstance(a, (Number,
                                                NumberSymbol,
                                                ImaginaryUnit)),
                        atom[:]))) \
        or (isinstance(atom, Pow) and isinstance(atom.base, (Number,
                                                             NumberSymbol,
                                                             ImaginaryUnit))
            and isinstance(atom.exp, (Number, NumberSymbol, ImaginaryUnit))):
        return [atom]
    elif isinstance(atom, (Add, Mul)):
        for a in atom:
            result.append(get_numbers(a))
    return result

def integer_divisors(n):
    n = abs(n)
    r = []
    for i in range(1, n/2+1):
        if n % i == 0:
            r.append(i)
    r.append(n)
    return r
