"""
translation of rmpoly , using monomials instead of packed exponents
for the monomials. Monomials are accessed via functions
with monomial_ prefix; their implementation can be a tuple
(or an ETuple ?)
"""

from sympy.polys.monomialtools import (
    monomial_mul, monomial_div, monomial_min,
)

from sympy.core import Add, Mul, Pow, Rational
from sympy.functions import sin, cos, exp, log
from sympy.polys.domains import QQ

from copy import copy

import re
import math

class TaylorEvalError(TypeError):
    pass

def giant_steps(target):
    """
    list of precision steps for the Newton's method

    code adapted from mpmath/libmp/libintmath.py
    """
    L = [target]
    while L[-1] > 2:
        L = L + [(L[-1]+1)//2]
    return L[::-1]

rpm = re.compile('[+-]')

def monomial_zero(n):
    return (0, )*n

def monomial_basis(i, n):
    """
    return the ith-basis element
    """
    a = [0]*n
    a[i] = 1
    return tuple(a)


def monomial_from_sequence(a):
    """
    return a monomial tuple from a list or tuple
    """
    return tuple(a)

def monomial_pow(a, n):
    """
    return the n-th pow of the monomial
    """
    b = [x*n for x in a]
    return tuple(b)

def monomial_tobasic(monom, *gens):
    term = []
    for g, m in zip(gens, monom):
        term.append(Pow(g, m))
    return Mul(*term)

class BaseLPoly(object):
    """
    abstract base class for polynomial rings of a ring K
    subclasses:
    RPoly polynomial ring on a commutative ring
          with no elements with more than one term
    MRPoly polynomial ring on a commutative ring
          with elements with more than one term
    NCPoly polynomial ring with noncommutative base ring
    The objects lp in one of these polynomial rings
    construct elements of Poly, a multivariate polynomial ring with
    monomial ordering.
    lp.K generates the elements in K
    K(0)  zero element in K
    K(1)  unit element in K
    O order object
    The string representation of polynomials is in decreasing order;
    in it ^ instead of ** is used as exponent symbol.

    >>> from sympy.core.numbers import Rational as QQ
    >>> from sympy.polys.monomialtools import lex
    >>> from sympy.polys.lpoly import LPoly

    >>> lp = LPoly('x, y', QQ, lex)
    >>> lp.zero_mon
    (0, 0)
    >>> x, y = lp.gens()
    >>> p = (x + y)**2
    >>> str(p)
     ' +x^2 +2*x*y +y^2'
    """

    def __init__(self, pol_gens, ring, O, **kwds):
        if isinstance(pol_gens, str):
            pol_gens = pol_gens.split(', ')
        self.pol_gens = tuple(pol_gens)
        self.ngens = len(pol_gens)
        self.ring = ring
        self.order = O
        if 'base_ring' in kwds:
            self.base_ring = kwds['base_ring']
        else:
            self.base_ring = None
        self.gens_dict = {}
        for i in range(self.ngens):
            self.gens_dict[self.pol_gens[i]] = i
        if 'parens' in kwds:
            self.parens = kwds['parens']
        else:
            self.parens = False
        if 'commuting' in kwds:
            self.commuting = kwds['commuting']
        else:
            self.commuting = True
        if isinstance(ring, BaseLPoly):
            self.parens = True
            if not ring.commuting:
                self.commuting = False
        str_ring = str(self.ring)
        if 'mpc' in str_ring or 'Complex' in str_ring:
            self.parens = True
        self.SR = False
        if hasattr(ring, '__name__'):
            if ring.__name__ == 'sympify':
                self.SR = True
                self.parens = True
        self.zero_mon = monomial_zero(self.ngens)

    def __str__(self):
        try:
            ring_name = self.ring.__name__
        except:
            ring_name = str(self.ring)
        s = 'LPoly with ngens=%d ring=%s' % (self.ngens, ring_name)

        return s

    def gens(self):
        """return the list of the variables
        """
        one = self.ring(1)
        ngens = self.ngens
        a = []
        for i in range(ngens):
            expv = monomial_basis(i, ngens)
            p = Poly(self)
            p[expv] = one
            a.append(p)
        return a

    def read_monom(self, s):
        """compute expv for the string 'x_0^e_0*...'
        where expv is (e_0, ..)
        """
        gens_dict = self.gens_dict
        ngens = self.ngens
        ring = self.ring
        s = s.strip()
        a = s.split('*')
        expv = [0]*ngens
        for x in a:
            t = x.split('^')
            ind = gens_dict[t[0]]
            if len(t) == 2:
                pw = int(t[1])
            else:
                pw = 1
            expv[ind] = pw
        return monomial_from_sequence(expv)

    def mon_eval(self, s):
        """compute the tuple (expv, c) for a string 'c*x_0^e_0*...'
        where expv = (e0, ..)
        """
        gens_dict = self.gens_dict
        ngens = self.ngens
        ring = self.ring
        s = s.strip()
        if s[0].isdigit() and '*' not in s:
            return (self.zero_mon, ring(s))
        if s[0].isdigit():
            a = s.split('*', 1)
            coeff = ring(a[0])
            a = a[1].split('*')
        else:
            a = s.split('*')
            coeff = ring(1)
        expv = [0]*ngens
        # e.g. a = ['x0^6', 'x1^6', 'x10^2']
        # split each element in ind, pw
        for x in a:
            t = x.split('^')
            ind = gens_dict[t[0]]
            if len(t) == 2:
                pw = int(t[1])
            else:
                pw = 1
            expv[ind] = pw
        return (tuple(expv), coeff)

    def read(self, s):
        """
        return the polynomial of s, in case of polynomial over polynomial
        s string of the form +(coeff1)*m1 + ... +(coeff0)
        where coeff0, coeff1, .. are polynomials in the field self.field
        m1, ..  monomial formed by the variables rp.gens()
        The form of the string must be exactly as above; each coeff is in brackets,
        there is always a + in front of the left bracket, if there is a monomial
        multiplying the monomial the right bracket is followed by * (no space allowed);
        the only exception is when s is a number.
        """
        s = s.strip()
        p = Poly(self)
        lp = p.lp
        zm = lp.zero_mon
        lp1 = self.ring
        if '(' not in s:
            p = Poly(p.lp)
            cp = lp1(s)
            if cp:
                p[zm] = cp
            return p
        assert s[0] == '+' and s[1] == '('
        # nparens = number of '(' parenthesis - number of ')' parenthesis
        nparens = 1
        # prev_pos starting position of the subexpression; the subexpression ends
        # when nparens becomes zero.
        prev_pos=2
        ns = len(s)
        pos = 2
        while 1:
            if s[pos] == '(':
                nparens += 1
            elif s[pos] == ')':
                nparens -= 1
            if nparens == 0:
                # found subexpression
                s1 = s[prev_pos:pos]
                # evaluate the subexpression
                cp = rp1(s1)
                if cp:
                    pos += 1
                    # if s1 is followed by *, then the monomial multiplying it follows
                    if pos < ns and s[pos] == '*':
                        pos += 1
                        # pos1= position of next coefficient subexpression
                        pos1 = s.find('+(', pos)
                        # +(
                        prev_pos = pos1+2
                        if pos1 == -1:
                            expv = self.read_monom(s[pos:])
                        else:
                            expv = self.read_monom(s[pos:pos1])
                    else:
                        expv = zm
                        pos1 = -1
                    p[expv] = cp
                    if pos1 == -1 or pos >= ns:
                        break
                    else:
                        pos = pos1 + 1
                    if pos >= ns:
                        break
                else:
                    # the coefficient is zero, go to the next term
                    if pos < ns:
                        pos1 = s.find('+(', pos)
                        if pos1 == -1:
                            break
                        prev_pos = pos1+2
                        pos = prev_pos
            else:
                pos += 1
                if pos >= ns:
                    break

        if nparens != 0:
            raise ValueError('nparens=%d' % nparens)
        return p



    def __call__(self, s):
        """
        generate a polynomial from a string; the string must be in
        the form c_0*m_0 + c_1*m_1 + ...
        where c_i are coefficients and m_i are monomials

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = lp('1/2 + 2/3*x^4 + 1/5*x*y^2')
        >>> str(p)
        ' +2/3*x^4 +1/5*x*y^2 +1/2'
        """
        if isinstance(s, Poly):
            if s.lp == self:
                return s
            else:
                raise NotImplementedError('cannot convert polynomial')
        if not isinstance(s, str):
            try:
                c = self.ring(s)
                p = Poly(self)
                if c:
                    p[self.zero_mon] = c
                return p
            except:
                raise ValueError('s cannot be converted to self.ring')
        if self.parens:
            return self.read(s)
        s = s.strip()
        a = []
        if (s[0] == '+' or s[0] == '-'):
            sgn = s[0]
        else:
            sgn = '+'
            s = '+' + s
        it = rpm.finditer(s)
        next(it)
        s1 = s[1:]
        if rpm.search(s1):
            start = 1
            for match in it:
                t=match.span()
                a.append((sgn, s[start:t[0]]))
                start=t[0]+1
                sgn = s[t[0]]
            a.append((sgn, s[t[1]:]))
        else:
            a.append((sgn, s1))
        p = Poly(self)
        for sgn, s1 in a:
            m = self.mon_eval(s1)
            if sgn == '-':
                m = (m[0], -m[1])
            if m[0] in p:
                p[m[0]] += m[1]
            else:
                p[m[0]] = m[1]
        p.strip()
        return p

    def from_mon(self, a):
        """polynomial from the monomial coeff*x_i^j
           a = (i, j, coeff)
        """
        p = Poly(self)
        ngens = self.ngens
        i = a[0]
        if i.__class__ == str:
            i = self.gens_dict[i]
        j = a[1]
        c = a[2]
        expv = [0]*ngens
        expv[i] = j
        expv = monomial_from_sequence(expv)
        p[expv] = c
        return p

    def from_expv(self, expv):
        "monomial from monomial tuple expv"
        p = Poly(self)
        p[expv] = 1
        return p




class LPoly(BaseLPoly):
    def __init__(self, pol_gens, ring, order, **kwds):
        BaseLPoly.__init__(self, pol_gens, ring, order, **kwds)
        self.__name__ = 'LPoly'


class MLPoly(BaseLPoly):
    """class of polynomials on a ring with elements with more than one term

    This class differs from LPoly only because the coefficients of
    its polynomials are surrounded by parenthesis

    >>> from sympy.polys.monomialtools import lex
    >>> from sympy.polys.lpoly import lgens, MLPoly
    >>> from sympy.core import sympify
    >>> from sympy.functions.elementary.exponential import log
    >>> lp, x, y = lgens('x, y', sympify, lex)
    >>> p = (x + x*log(2) + y)**2
    >>> lp.__class__
     <class 'sympy.polys.lpoly.MLPoly'>
    """
    def __init__(self, pol_gens, ring, order, **kwds):
        BaseRPoly.__init__(self, pol_gens, ring, order, **kwds)
        self.parens = True
        self.__name__ = 'MLPoly'


class NCLPoly(BaseLPoly):
    """class of polynomials on noncommuting ring
    """
    def __init__(self, pol_gens, ring, order, **kwds):
        BaseRPoly.__init__(self, pol_gens, ring, order, **kwds)
        self.parens = True
        self.commuting = False
        self.__name__ = 'NCLPoly'

def lgens(pol_gens, ring, order, **kwds):
    """
     factory function to generate LPoly object and its generators

    >>> from sympy.core.numbers import Rational as QQ
    >>> from sympy.polys.monomialtools import lex
    >>> from sympy.polys.lpoly import lgens
    >>> lp, x, y = lgens('x, y', QQ, lex)
    >>> p = (x + y)**2
    >>> str(p)
    ' +x^2 +2*x*y +y^2'
    """
    lp = BaseLPoly(pol_gens, ring, order, **kwds)
    if lp.parens:
        lp.__class__ = MLPoly
        lp.__name__ = 'MLPoly'
    else:
        lp.__class__ = LPoly
        lp.__name__ = 'LPoly'
    return [lp] + lp.gens()

def mlgens(pol_gens, ring, **kwds):
    """
    factory function to generate MLPoly object and its generators
    """
    lp = MLPoly(pol_gens, ring, **kwds)
    return [lp] + lp.gens()

def ncrgens(pol_gens, bits_exp, ring, **kwds):
    """
    factory function to generate NCRPoly object and its generators
    """
    lp = NCLPoly(pol_gens, ring, **kwds)
    return [lp] + lp.gens()




class Poly(dict):
    """
    elements of a multivariate polynomial ring

    >>> from sympy.core.numbers import Rational as QQ
    >>> from sympy.polys.monomialtools import lex
    >>> from sympy.polys.lpoly import lgens, Poly
    >>> lp, x, y = lgens('x, y', QQ, lex)
    >>> p = Poly(lp)
    >>> str(p)
    '0'
    >>> p1 = x + y**2
    >>> isinstance(p1, Poly)
    True
    """
    def __init__(self, lp, **kw):
        self.lp = lp
        dict.__init__(self, **kw)

    def __str__(self):
        """string representation of a polynomial
        Terms are in decreasing order; term coefficients are in front

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = x**2/3 + 2*x*y**2/5
        >>> str(p)
        ' +1/3*x^2 +2/5*x*y^2'
        """
        if not self:
            return '0'
        lp = self.lp
        if lp.parens:
            return self.tostr()
        pol_gens = lp.pol_gens
        ngens = lp.ngens
        zm = self.lp.zero_mon
        s = ''
        a = list(self.keys())
        a.sort(key=lambda m: lp.order(m), reverse=True)
        z = monomial_zero(ngens)
        for expv in a:
            c = self[expv]
            if c > 0:
                s += ' +'
            else:
                s += ' -'
            if c < 0:
                c = -c
            if c != 1 or expv == zm:
                cnt1 = str(c)
            else:
                cnt1 = ''
            sa = []
            for i in range(ngens):
                exp = expv[i]
                if exp > 1:
                    sa.append('%s^%d' % (pol_gens[i], exp))
                if exp == 1:
                    sa.append('%s' % pol_gens[i])
            if cnt1:
                sa = [cnt1] + sa
            s += '*'.join(sa)
        return s

    __repr__ = __str__

    def tostr(self):
        lp = self.lp
        pol_gens = lp.pol_gens
        ngens = lp.ngens
        zm = lp.zero_mon
        s = ''
        a = list(self.keys())
        a.sort(key=lambda m: lp.order(m), reverse=True)
        for expv in a:
            c = self[expv]
            s += ' +(%s)' % c
            if expv != zm:
                s += '*'
            # throw away the bits for the hidden variable
            i = 0
            sa = []
            for i in range(ngens):
                exp = expv[i]
                if exp > 1:
                    sa.append('%s^%d' % (pol_gens[i], exp))
                if exp == 1:
                    sa.append('%s' % pol_gens[i])
            s += '*'.join(sa)

        return s

    def strip(p):
        """eliminate monomials with zero coefficient
        """
        z = p.lp.ring(0)
        for k, v in list(p.items()):
            if v == z:
                del p[k]

    def variables(self):
        """return the tuple of the indices of the
        variables occurring in the polynomial p

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y, z = lgens('x, y, z', QQ, lex)
        >>> p = x + y**2 + x*y
        >>> p.variables()
        (0, 1)
        """
        a = []
        for expv in self:
            if(self[expv] == 0):
                continue
            i = 0
            for i, e in enumerate(expv):
                if e and i not in a:
                    a.append(i)
        a.sort()
        return tuple(a)

    def __eq__(p1, p2):
        """
        equality test for polynomials

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p1 = (x+y)**2 + (x-y)**2
        >>> p2 = 4*x*y
        >>> p1 == p2
        False
        >>> p2 = 2*(x**2 + y**2)
        >>> p1 == p2
        True
        """
        if not p2 or p2 == 0:
            return not p1
        lp1 = p1.lp
        if isinstance(p2, Poly):
            if lp1 == p2.lp:
                return dict.__eq__(p1, p2)
        else:
            zm = lp1.zero_mon
            # assume p2 is a coefficient
            if zm not in p1 or len(p1) > 1:
                return False
            return p1[zm] == lp1.ring(p2)

    def __add__(p1, p2):
        """add two polynomials

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = (x+y)**2 + (x-y)**2
        >>> str(p)
        ' +2*x^2 +2*y^2'
        """
        if not p2:
            return p1.copy()
        lp1 = p1.lp
        zm = lp1.zero_mon
        if isinstance(p2, Poly):
            if lp1 == p2.lp:
                p = Poly(p1.lp)
                for k, v in p1.items():
                    if k in p2:
                        r = v + p2[k]
                        if r:
                            p[k] = r
                    else:
                        p[k] = v
                for k, v in p2.items():
                    if k not in p1:
                        p[k] = v
                return p
            elif p2.lp == lp1.ring:
                p = p1.copy()
                if zm not in list(p1.keys()):
                    p[zm] = lp1.ring(p2)
                else:
                    if p2 == -p[zm]:
                        del p[zm]
                    else:
                        p[zm] += p2
                return p
            elif lp1 == p2.lp.ring:
                return p2+p1
            else:
                raise ValueError('cannot sum p1 and p2')
        # assume p2 in a number
        else:
            p = p1.copy()
            cp2 = lp1.ring(p2)
            if not cp2:
                return p
            if zm not in list(p1.keys()):
                p[zm] = cp2
            else:
                if p2 == -p[zm]:
                    del p[zm]
                else:
                    p[zm] += cp2
            return p

    def copy(self):
        """
        copy of a polynomial
        polynomials are mutable; if one is interested in preserving
        a polynomial, and one plans to use inplace operations, one
        can copy the polynomial

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = (x+y)**2
        >>> p1 = p.copy()
        >>> p2 = p
        >>> p[lp.zero_mon] = 3
        >>> str(p)
        ' +x^2 +2*x*y +y^2 +3'
        >>> str(p1)
        ' +x^2 +2*x*y +y^2'
        >>> str(p2)
        ' +x^2 +2*x*y +y^2 +3'
        """
        return copy(self)

    def coefficient(self, p1):
        """the coefficient of a monomial p1 is the sum ot the terms in
        self which have the same degrees in the variables present in p1,
        divided by p1

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = (1 + x + y)**3
        >>> p1 = p.coefficient(y**2)
        >>> str(p1)
        ' +3*x +3'
        """
        lp = self.lp
        order = lp.order
        k = list(p1.keys())
        if len(k) != 1:
            raise TypeError('the argument of coeff must be a monomial')
        expv1 = k[0]
        # mask1 used to select the exponents of the variables present in p1
        v1 = p1.variables()
        p = Poly(lp)
        zm = lp.zero_mon
        for expv in self:
            b = 1
            for i in v1:
                if expv[i] != expv1[i]:
                    b = 0
                    break
            if b:
                p[monomial_div(expv, expv1)] = self[expv]
        return p

    def coeff(self, p1):
        """monomial coefficient in self of the leading monomial of p1

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = (1 + x + y)**3/5
        >>> p.coeff(y**2)
        3/5
        >>> p.coeff(x*y)
        6/5
        """
        expv = p1.leading_expv()
        if expv in self:
            return self[expv]
        else:
            return lp.ring(0)

    def subs(p, **rules):
        """
        substitution
        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = x**2 + y**2
        >>> p1 = p.subs(x=x+y, y=x-y)
        >>> str(p1)
        ' +2*x^2 +2*y^2'
        """
        lp = p.lp
        sb = LPolySubs(lp, lp, rules)
        return sb.subs(p)

    def subs_trunc(p, iv, nv, **rules):
        """
        substitution with truncation

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = x**2 + y**2
        >>> p2 = p.subs_trunc('x', 3, y=(x+y)**2)
        >>> str(p2)
        ' +6*x^2*y^2 +x^2 +4*x*y^3 +y^4'
        """
        lp = p.lp
        sb = LPolySubs(lp, lp, rules)
        return sb.subs_trunc(p, iv, nv)


    def __radd__(p1, n):
        # assume n is in p1.lp.ring
        p = p1.copy()
        if not n:
            return p
        lp = p1.lp
        zm = lp.zero_mon
        if zm not in list(p1.keys()):
            p[zm] = lp.ring(n)
        else:
            if n == -p[zm]:
                del p[zm]
            else:
                p[zm] += n
        return p

    def __neg__(self):
        return self*(-1)

    def __pos__(self):
        return self

    def __iadd__(p1, p2):
        """inplace add polynomials

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p1 = x + y**2
        >>> p2 = x**2
        >>> p3 = p1
        >>> p1 += p2
        >>> str(p1)
        ' +x^2 +x +y^2'
        >>> p1 == p3
        True
        """
        if not p2:
            return p1
        lp1 = p1.lp
        if isinstance(p2, Poly):
            if lp1 != p2.lp:
                raise ValueError('p1 and p2 must have the same lp')
            dl = []
            for k, v in p1.items():
                if k in p2:
                    if p2[k] == -v:
                        dl.append(k)
                    else:
                        p1[k] = v + p2[k]
            for k in p2:
                if k not in p1:
                    p1[k] = p2[k]
            for k in dl:
                del p1[k]
            return p1
        else:
            mz = lp1.zero_mon
            if mz not in list(p1.keys()):
                p1[mz] = lp1.ring(p2)
            else:
                if p2 == -p1[mz]:
                    del p1[mz]
                else:
                    p1[mz] += p2
            return p1

    def __sub__(p1, p2):
        """
        subtraction of polynomials

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p1 = x + y**2
        >>> p2 = x*y + y**2
        >>> p3 = p1 - p2
        >>> str(p3)
        ' -x*y +x'
        """
        if not p2:
            return p1.copy()
        lp1 = p1.lp
        mz = lp1.zero_mon
        p = Poly(lp1)
        if isinstance(p2, Poly):
            if lp1 == p2.lp:
                for k in p1:
                    if k in p2:
                        r = p1[k] - p2[k]
                        if r:
                            p[k] = r
                    else:
                        p[k] = p1[k]
                for k in p2:
                    if k not in p1:
                        p[k] = -p2[k]
                return p
            elif p2.lp == lp1.ring:
                p = p1.copy()
                if mz not in list(p1.keys()):
                    p[mz] = -lp1.ring(p2)
                else:
                    if p2 == p[mz]:
                        del p[mz]
                    else:
                        p[mz] -= p2
                return p
            else:
                raise ValueError('cannot coerce p2')
        # assume p2 in a number
        else:
            p2 = lp1.ring(p2)
            p = copy(p1)
            if mz not in list(p1.keys()):
                p[mz] = -p2
            else:
                if p2 == p[mz]:
                    del p[mz]
                else:
                    p[mz] -= p2
            return p

    def __rsub__(p1, n):
        """
        n - p1 with n convertible to the coefficient ring

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = x + y
        >>> p1 = 4 - p
        >>> str(p1)
        ' -x -y +4'
        """
        p = Poly(p1.lp)
        for expv in p1:
            p[expv] = -p1[expv]
        p += n
        return p


    def __isub__(p1, p2):
        """
        inplace subtraction of polynomials
        """
        lp1 = p1.lp
        if isinstance(p2, Poly):
            if lp1 != p2.lp:
                raise ValueError('p1 and p2 must have the same lp')
            dl = []
            for k in p1:
                if k in p2:
                    if p2[k] == p1[k]:
                        dl.append(k)
                    else:
                        p1[k] = p1[k] - p2[k]
            for k in p2:
                if k not in p1:
                    p1[k] = -p2[k]
            for k in dl:
                del p1[k]
            return p1
        else:
            mz = lp1.zero_mon
            if mz not in p1:
                p1[mz] = -p2
            else:
                if p1[mz] == p2:
                    del p1[mz]
                else:
                    p1[mz] -= p2
            return p1

    def __mul__(p1, p2):
        """multiply two polynomials

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p1 = x + y
        >>> p2 = x - y
        >>> p3 = p1*p2
        >>> str(p3)
        ' +x^2 -y^2'
        """
        lp1 = p1.lp
        p = Poly(lp1)
        if not p2:
            return p
        if isinstance(p2, Poly):
            if lp1 == p2.lp:
                get = p.get
                p2it = list(p2.items())
                for exp1, v1 in p1.items():
                    for exp2, v2 in p2it:
                        exp = monomial_mul(exp1, exp2)
                        p[exp] = get(exp, 0) + v1*v2
                p.strip()
                return p
            if p2.lp != lp1.ring:
                if lp1 == p2.lp.ring:
                    p = Poly(p2.lp)
                    for exp2, v2 in p2.items():
                        p[exp2] = p1*v2
                    return p
                else:
                    raise ValueError('p1 and p2 must have the same lp')
        # assume p2 in a number
        for exp1, v1 in p1.items():
            v = v1*p2
            if v:
                p[exp1] = v
        return p

    def __rmul__(p1, p2):
        """
         p2 + p1 with p2 in the coefficient ring of p1

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = x + y
        >>> p1 = 4 * p
        >>> str(p1)
        ' +4*x +4*y'
        """
        p = Poly(p1.lp)
        if not isinstance(p2, Poly):
            if not p2:
                return p
        for exp1, v1 in p1.items():
            v = p2*v1
            if v:
                p[exp1] = v
        return p


    def mul_iadd(p, p1, p2):
        """p += p1*p2

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = x**2
        >>> p1 = x + y
        >>> p2 = x - y
        >>> p.mul_iadd(p1, p2)
        >>> str(p)
        ' +2*x^2 -y^2'
        """
        if isinstance(p1, Poly) and isinstance(p2, Poly):
            if p1.lp != p2.lp:
                raise ValueError('p1 and p2 must have the same lp')
            get = p.get
            for exp1, v1 in p1.items():
                for exp2, v2 in p2.items():
                    exp = monomial_mul(exp1, exp2)
                    p[exp] = get(exp, 0) + v1*v2
            p.strip()
        else:
            raise NotImplementedError

    def imul_num(p, c):
        """multiply inplace the polynomial p by an element in the
        coefficient ring

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = x + y**2
        >>> p.imul_num(3)
        >>> str(p)
        ' +3*x +3*y^2'
        """
        if not c:
            p.clear()
            return
        for exp in p:
            p[exp] *= c

    def mul_num(p, c):
        """multiply the polynomial p by an element in the
        coefficient ring

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = x + y**2
        >>> p1 = p.mul_num(3)
        >>> str(p1)
        ' +3*x +3*y^2'
        """
        if not c:
            return p.lp(0)
        p1 = p.copy()
        for exp in p:
            p1[exp] = p[exp]*c
        return p1

    def __truediv__(p1, p2):
        """division by a term in the coefficient ring or
        exact division by a polynomial

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p1 = (x**2 + x + y)*(x**2-y**2)
        >>> p2 = x + y
        >>> p3 = p1/p2
        >>> p4 = (x**2 + x + y)*(x-y)
        >>> p3 == p4
        True
        """
        lp1 = p1.lp
        if isinstance(p2, Poly):
            if lp1 == p2.lp:
                q, r = p1.division([p2])
                if r:
                    raise NotImplementedError('__div__ performs only division without remainder')
                return q[0]
            elif p2.lp == lp1.ring:
                zm = p2.lp.zero_mon
                p = Poly(lp1)
                # if p is not a constant, not implemented
                if list(p2.keys()) != [zm]:
                    raise NotImplementedError
                else:
                    p2 = p2[zm]
            else:
                raise NotImplementedError('cannot divide p1 by p2')
        # assume p2 in a number
        p = Poly(lp1)
        if not p2:
            raise ZeroDivisionError
        try:
            for exp1, v in p1.items():
                p[exp1] = v/p2
        except:
            if not lp1.base_ring:
                raise NotImplementedError('cannot divide p1 by p2')
            else:
                c = lp1.base_ring(1)/p2
                for exp1, v in p1.items():
                    p[exp1] = v*c
        return p

    __div__ = __truediv__

    def iadd_mon(self, a):
        """add inplace the monomial coeff*x0^i0*x1^i1*...
        a = ((i0, i1, ...), coeff)

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> from sympy.polys.lpoly import monomial_from_sequence
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = x**4 + 2*y
        >>> m = monomial_from_sequence((1, 2))
        >>> p.iadd_mon((m, 5))
        >>> str(p)
        ' +x^4 +5*x*y^2 +2*y'
        """
        coeff = a[1]
        expv = a[0]
        if expv in self:
            self[expv] += coeff
            if self[expv] == 0:
                del self[expv]
        else:
            self[expv] = coeff

    def iadd_m_mul_q(p1, p2, xxx_todo_changeme):
        """ p1 += p2*monom
        m monomial tuple
        c coefficient

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> from sympy.polys.lpoly import monomial_from_sequence
        >>> lp, x, y, z = lgens('x, y, z', QQ, lex)
        >>> p1 = x**4 + 2*y
        >>> p2 = y + z
        >>> m = monomial_from_sequence((1, 2, 3))
        >>> p3 = p1.iadd_m_mul_q(p2, (m, 3))
        >>> str(p1)
        ' +x^4 +3*x*y^3*z^3 +3*x*y^2*z^4 +2*y'
        >>> p1 == p3
        True
        """
        (m, c) = xxx_todo_changeme
        get = p1.get
        c = p1.lp.ring(c)
        for k, v in p2.items():
            ka = monomial_mul(k, m)
            p1[ka] = get(ka, 0) + v*c
        p1.strip()
        return p1

    def leading_expv(self):
        """leading monomial tuple according to the monomial ordering

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y, z = lgens('x, y, z', QQ, lex)
        >>> p = x**4 + x**3*y + x**2*z**2 + z**7
        >>> p.leading_expv()
        (4, 0, 0)
        """
        if self:
            order = self.lp.order
            return max(self, key=lambda m: order(m))
        else:
            return None

    def leading_term(self):
        """leading term according to the monomial ordering

       >>> from sympy.core.numbers import Rational as QQ
       >>> from sympy.polys.monomialtools import lex
       >>> from sympy.polys.lpoly import lgens
       >>> lp, x, y = lgens('x, y', QQ, lex)
       >>> p = (x+y)**4
       >>> p1 = p.leading_term()
       >>> str(p1)
       ' +x^4'
       """
        p = Poly(self.lp)
        expv = self.leading_expv()
        if expv:
            p[expv] = self[expv]
        return p

    def square(p1):
        """square of a polynomial

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = x + y**2
        >>> p1 = p.square()
        >>> str(p1)
        ' +x^2 +2*x*y^2 +y^4'
        """
        lp = p1.lp
        if not lp.commuting:
            return p1*p1
        p = Poly(lp)
        get = p.get
        keys = list(p1.keys())
        for i in range(len(keys)):
            k1 = keys[i]
            pk = p1[k1]
            for j in range(i):
                k2 = keys[j]
                exp = monomial_mul(k1, k2)
                p[exp] = get(exp, 0) + pk*p1[k2]
        p.imul_num(2)
        for k, v in p1.items():
            k2 = monomial_mul(k, k)
            p[k2] = get(k2, 0) + v**2
        p.strip()
        return p

    def pow_miller(p, m):
        """power of an univariate polynomial

        p = sum_i=0^L p_i*x**i
        p**m = sum_k=0^(m*L) a(m,k)*x**k
        Miller pure recurrence formula (see article by
        D. Zeilberger)
        a(m,k) = 1/(k*p_0)*sum_i=1^L p_i*((m+1)*i-k)*a(m,k-i)

        Reference: D. Zeilberger 'The Miller Recurrence for
        Exponentatiating a Polynomial, and its q-Analog'
        """

        lp = p.lp
        degp = max(p)[0]
        pv = [0]*(degp+1)
        for k,v in p.iteritems():
            pv[k[0]] = v
        a = [lp(0) for i in range(m*degp+1)]
        a[0] = pv[0]**m
        res = lp(a[0])
        for k in range(1, m*degp+1):
            s = 0
            for i in range(1, min(degp,k)+1):
                s += pv[i]*((m+1)*i - k)*a[k-i]
            a[k] = s/(k*pv[0])
            res[(k,)] = a[k]
        return res

    def __pow__(self, n):
        """power of a polynomial

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = x + y**2
        >>> p1 = p**3
        >>> str(p1)
        ' +x^3 +3*x^2*y^2 +3*x*y^4 +y^6'
        """
        lp = self.lp
        # test if n is an integer
        if not isinstance(n, int):
            raise NotImplementedError
        n = int(n)
        if n < 0:
            raise ValueError('n >= 0 is required')
        if n == 0:
            if self:
                return lp(1)
            else:
                raise ValueError
        elif len(self) == 1:
            p = Poly(lp)
            k, v = list(self.items())[0]
            # treat case abs(v) = 1 separately to deal with the case
            # in which n is too large to be allowed in v**n
            kn = monomial_pow(k, n)
            if v == 1:
                p[kn] = v
            elif v == -1:
                if n%2 == 0:
                    p[kn] = -v
                else:
                    p[kn] = v
            else:
                p[kn] = v**n
            return p
        elif n == 1:
            return copy(self)
        elif n == 2:
            return self.square()
        elif n == 3:
            return self*self.square()
        p = lp(1)
        while 1:
            if n&1:
                p = p*self
                n -= 1
                if not n:
                    break
            self = self.square()
            n = n // 2
        return p

    def division(self, fv):
        """division algorithm, see [CLO] p64
        fv array of polynomials
           return qv, r such that
           self = sum(fv[i]*qv[i]) + r

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> f = x**3
        >>> f0 = x-y**2
        >>> f1 = x-y
        >>> qv, r = f.division((f0, f1))
        >>> str(qv[0])
        ' +x^2 +x*y^2 +y^4'
        >>> str(qv[1])
        '0'
        >>> str(r)
        ' +y^6'
        """
        lp = self.lp
        if not self:
            return [], Poly(lp)
        for f in fv:
            if f.lp != lp:
                raise ValueError('self and f must have the same lp')
        s = len(fv)
        qv = [Poly(lp) for i in range(s)]
        p = self.copy()
        r = Poly(lp)
        order = self.lp.order
        expvs = [fx.leading_expv() for fx in fv]
        while p:
            i = 0
            divoccurred = 0
            while i < s and divoccurred == 0:
                expv = p.leading_expv()
                expv1 = monomial_div(expv, expvs[i])
                if expv1:
                    c = p[expv]/fv[i][expvs[i]]
                    qv[i].iadd_mon((expv1, c))
                    p.iadd_m_mul_q(fv[i], (expv1, -c))
                    divoccurred = 1
                else:
                    i += 1
            if not divoccurred:
                expv =  p.leading_expv()
                r.iadd_mon((expv, p[expv]))
                del p[expv]
        if expv == lp.zero_mon:
            r += p
        return qv, r

    def mod1(self, fv):
        """fv array of (expv, p), p monic
           return r such that
           self = sum(fv[i]*qv[i]) + r
        """
        lp = self.lp
        if not self:
            return Poly(lp)
        s = len(fv)
        p = self.copy()
        r = Poly(lp)
        expvs = [fx[0] for fx in fv]
        fv = [fx[1] for fx in fv]
        expv = p.leading_expv()
        while p:
            i = 0
            divoccurred = 0
            while i < s and divoccurred == 0:
                expv1 = monomial_div(expv, expvs[i])
                if expv1:
                    p.iadd_m_mul_q(fv[i], (expv1, -p[expv]))
                    expv = p.leading_expv()
                    divoccurred = 1
                else:
                    i += 1
            if not divoccurred:
                r.iadd_mon((expv, p[expv]))
                del p[expv]
                expv = p.leading_expv()
        if expv == lp.zero_mon:
            r += p
        return r

    def trunc(p1, i, h):
        """monomials containing x^k, k >= h neglected
        i is the name of the variable x, or its index

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = (x+y)**4
        >>> p1 = p.trunc('x', 3)
        >>> str(p1)
        ' +6*x^2*y^2 +4*x*y^3 +y^4'
        """
        lp = p1.lp
        if isinstance(h,int):
            h = int(h)
            if isinstance(i, str):
                i = lp.gens_dict[i]
            p = Poly(lp)
            for exp1 in p1:
                if exp1[i] >= h:
                    continue
                p[exp1] = p1[exp1]
            return p
        elif isinstance(h, tuple) or isinstance(h, list):
            raise NotImplementedError


    def mul_trunc(p1, p2, i, h):
        """truncation of p1*p2
        p1 and p2 polynomials
        If h is an integer,
        let i is the name of the variable x, or its index;
        neglect in p1*p2 the monomials containing x^k, k >= h
        If h is a tuple or list of integers,
        i_j is the name of the variable x_j, or its index, for each
        i_j in i; h_j is the corresponding power in h
        neglect in p1*p2 the monomials containing x_j^k, k >= i_j

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p1 = (x + y)**2
        >>> p2 = p1.mul_trunc(x**2, 'x', 3)
        >>> str(p2)
        ' +x^2*y^2'
        """
        lp = p1.lp
        if p1.lp != p2.lp:
            raise ValueError('p1 and p2 must have the same lp')
        if isinstance(h, int):
            h = int(h)
            if isinstance(i, str):
                iv = lp.gens_dict[i]
            else:
                iv = i
            p = Poly(p1.lp)
            get = p.get
            items2 = list(p2.items())
            items2.sort(key=lambda e: e[0][iv])
            for exp1, v1 in p1.items():
                for exp2, v2 in items2:
                    if exp1[iv] + exp2[iv] < h:
                        exp = monomial_mul(exp1, exp2)
                        p[exp] = get(exp, 0) + v1*v2
                    else:
                        break
            p.strip()
            return p
        elif isinstance(h, tuple) or isinstance(h, list):
            raise NotImplementedError

    def square_trunc(p1, i, h):
        """truncation of p1*p1
        If h is an integer, let i is the name of the variable x, or its index;
        neglect in p1*p1 the monomials containing x^k, k >= h
        If h is a tuple or list of integers,
        i_j is the name of the variable x_j, or its index, for each
        i_j in i; h_j is the corresponding power in h
        neglect in p1*p2 the monomials containing x_j^k, k >= i_j

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p1 = (x + y)**2
        >>> p2 = p1.square_trunc('x', 3)
        >>> str(p2)
        ' +6*x^2*y^2 +4*x*y^3 +y^4'
        """
        lp = p1.lp
        if not lp.commuting:
            return p1.mul_trunc(p1, i, h)
        if isinstance(h, int):
            h = int(h)
            if isinstance(i, str):
                iv = lp.gens_dict[i]
            else:
                iv = i
            p = Poly(lp)
            get = p.get
            items = list(p1.items())
            items.sort(key=lambda e: e[0][iv])
            for i in range(len(items)):
                exp1, v1 = items[i]
                for j in range(i):
                    exp2, v2 = items[j]
                    if exp1[iv] + exp2[iv] < h:
                        exp = monomial_mul(exp1, exp2)
                        p[exp] = get(exp, 0) + v1*v2
                    else:
                        break
            p.imul_num(2)
            for expv, v in p1.items():
                if 2*expv[iv] < h:
                    e2 = monomial_mul(expv, expv)
                    p[e2] = get(e2, 0) + v**2
            p.strip()
            return p
        elif isinstance(h, tuple) or isinstance(h, list):
            raise NotImplementedError

    def pow_trunc(self, n, i, h):
        """truncation of self**n using mul_trunc

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p1 = 1 + x + y
        >>> p2 = p1.pow_trunc(3, 'y', 2)
        >>> str(p2)
        ' +x^3 +3*x^2*y +3*x^2 +6*x*y +3*x +3*y +1'
        """
        lp = self.lp
        if n != int(n):
            raise NotImplementedError
        n = int(n)
        if n == 0:
            if self:
                return lp(1)
            else:
                raise ValueError
        if n < 0:
            # for univariate series series_inversion is faster than pow)trunc
            if isinstance(h,int):
                p1 = self.pow_trunc(-n, i, h)
                return p1.series_inversion(i, h)
            else:
                p1 = self.series_inversion(i, h)
                return p1.pow_trunc(-n, i, h)
        if n == 1:
            return self.trunc(i, h)
        if n == 2:
            return self.square_trunc(i, h)
        if n == 3:
            p2 = self.square_trunc(i, h)
            return p2.mul_trunc(self, i, h)
        p = lp(1)
        while 1:
            if n&1:
                p = self.mul_trunc(p, i, h)
                n -= 1
                if not n:
                    break
            self = self.square_trunc(i, h)
            n = n // 2
        return p

    def has_constant_term(p, iv):
        lp = p.lp
        if isinstance(iv, str):
            iv = [lp.gens_dict[iv]]
        else:
            if isinstance(iv[0], str):
                d = lp.gens_dict
                iv = [d[x] for x in iv]
        zm = lp.zero_mon
        a = [0]*lp.ngens
        for i in iv:
            a[i] = 1
        miv = monomial_from_sequence(a)
        for expv in p:
            if monomial_min(expv, miv) == zm:
                return True
        return False

    def _series_inversion1(p, iv, nv):
        """univariate series inversion 1/p
        iv name of the series variable
        nv precision of the series

        The Newton method is used.
        """
        lp = p.lp
        zm = lp.zero_mon
        if zm not in p:
            raise ValueError('no constant term in series')
        if (p - p[zm]).has_constant_term(iv):
            raise ValueError('p cannot contain a constant term depending on parameters')
        if p[zm] != lp.ring(1):
            # TODO add check that it is a unit
            p1 = lp(1)/p[zm]
        else:
            p1 = lp(1)
        for prec in giant_steps(nv):
            tmp = p1.square()
            tmp = tmp.mul_trunc(p, iv, prec)
            p1 = 2*p1 - tmp
        return p1


    def series_inversion(p, iv, nv):
        """multivariate series inversion 1/p
        iv list of variable names or variable indices
        nv list of truncations for these variables
        In the case of one variable it can also be:
          iv variable name or variable index (0)
          nv truncation integer for the variable
        p is a series with O(x_1^n_1*..x_m^n_m) in
        variables x_k with index or name iv[k-1]
        p has constant term different from zero

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> from sympy.core import sympify
        >>> from sympy import Symbol
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = (1 + y*x).series_inversion('x', 6)
        >>> str(p)
        ' -x^5*y^5 +x^4*y^4 -x^3*y^3 +x^2*y^2 -x*y +1'
        >>> a = Symbol('a')
        >>> lp, x = lgens('x', sympify, lex)
        >>> p1 = (1 + a + x).series_inversion('x', 3)
        >>> str(p1)
        ' +((a + 1)**(-3))*x^2 +(-1/(a + 1)**2)*x +(1/(a + 1))'
        """
        lp = p.lp
        zm = lp.zero_mon
        if zm not in p:
            raise NotImplementedError('no constant term in series')
        if (p-p[zm]).has_constant_term(iv):
            raise NotImplementedError('p-p[0] must not have a constant term in the series variables')
        if not lp.commuting:
            return p._series_inversion_nc(iv, nv)
        if isinstance(nv, int):
            return p._series_inversion1(iv, nv)
        raise NotImplementedError

    def derivative(self, n):
        """derivative of p with respect to x_n; the argument n is the
        index of the variable x_n = self.lp.gens()[n], or the
        corresponding string.

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = x + x**2*y**3
        >>> p1 = p.derivative('x')
        >>> str(p1)
        ' +2*x*y^3 +1'
        """
        lp = self.lp
        pol_gens = lp.pol_gens
        if n.__class__ == str:
            n = lp.gens_dict[n]
        else:
            n = int(n)
        if n.__class__ == str:
            n = lp.gens_dict[n]
        p1 = Poly(lp)
        mn = monomial_basis(n, lp.ngens)
        for expv in self:
            if expv[n]:
                e = monomial_div(expv, mn)
                p1[e] = self[expv]*expv[n]
        return p1


    def integrate(self, n):
        """ integrate p with respect to x_n; the argument n is the
        index of the variable x_n = self.lp.gens()[n], or the
        corresponding string.

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = x + x**2*y**3
        >>> p1 = p.integrate('x')
        >>> str(p1)
        ' +1/3*x^3*y^3 +1/2*x^2'
        """
        lp = self.lp
        pol_gens = lp.pol_gens
        if n.__class__ == str:
            n = lp.gens_dict[n]
        p1 = Poly(lp)
        mn = monomial_basis(n, lp.ngens)
        for expv in self:
            e = monomial_mul(expv, mn)
            p1[e] = self[expv]/(expv[n]+1)
        return p1

############## elementary functions ####################

    def series_from_list(p, c, iv, nv, concur=1):
        """series sum c[n]*p^n
        reduce the number of multiplication summing concurrently
        ax = [1, p, p^2, .., p^(J-1)]
        s = sum(c[i]*ax[i] for i in range(0, J)) +
            sum(c[i]*ax[i] for i in range(J, 2*J))*p^J +
            sum(c[i]*ax[i] for i in range(2*J, 3*J))*p^(2*J) + ...+
            sum(c[i]*ax[i] for i in range((K-1)*J, K*J))*p^((K-1)*J)
        with K >= (n+1)/J

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> c = [1, 3, 5, 7]
        >>> p1 = (x+y).series_from_list(c, 'x', 3)
        >>> p2 = (1 + 3*(x+y) + 5*(x+y)**2 + 7*(x+y)**3).trunc('x', 3)
        >>> p1 == p2
        True
        """
        lp = p.lp
        zm = lp.zero_mon
        assert zm not in p
        n = len(c)
        if not concur:
            q = lp(1)
            s = c[0]*q
            for i in range(1, n):
                q = q.mul_trunc(p, iv, nv)
                s += c[i]*q
            return s
        J = int(math.sqrt(n) + 1)
        if n%J == 0:
            K = n//J
        else:
            K = int(n/J) + 1
        ax = [lp(1)]
        b = 1
        q = lp(1)
        # TODO get a good value
        if len(p) < 20:
            for i in range(1, J):
                q = q.mul_trunc(p, iv, nv)
                ax.append(q)
        else:
            for i in range(1, J):
                if i%2 == 0:
                    q = ax[i//2].square_trunc(iv, nv)
                else:
                    q = q.mul_trunc(p, iv, nv)
                ax.append(q)
        # optimize using square_trunc
        pj = ax[-1].mul_trunc(p, iv, nv)
        b = lp(1)
        s = lp(0)
        for k in range(K-1):
            r = J*k
            s1 = c[r]
            for j in range(1, J):
                s1 += c[r+j]*ax[j]
            s1 = s1.mul_trunc(b, iv, nv)
            s += s1
            b = b.mul_trunc(pj, iv, nv)
            if not b:
                break
        k = K-1
        r = J*k
        if r < n:
            s1 = c[r]*lp('1')
            for j in range(1, J):
                if r+j >= n:
                    break
                s1 += c[r+j]*ax[j]
            s1 = s1.mul_trunc(b, iv, nv)
            s += s1
        return s

    def fun(p, f, *args):
        """
        function of a multivariate series computed by substitution

        p multivariate series
        f method name or function
        args[:-2] arguments of f, apart from the first one
        args[-2] = iv names of the series variables
        args[-1] = nv list of the precisions of the series variables
        The case with f method name is used to compute tan and nth_root
        of a multivariate series:
        p.fun('tan', iv, nv)
        compute _x.tan(iv, sum(nv)), then substitute _x with p
        p.fun('nth_root', n, iv, nv)
        compute (_x + p[0]).nth_root(n, '_x', sum(nv)), then substitute _x
        with p - p[0]
        example with f function:
        f = Poly.exp_series0
        p1 = p.fun(Poly.exp_series0, ['y', 'x'], [h, h])
        """
        lp = p.lp
        lp1 = LPoly(['_x'], lp.ring, lp.order)
        _x = lp1.gens()[0]
        h = args[-1]
        if not isinstance(h,int):
            h = sum(args[-1])
        args1 = args[:-2] + ('_x', h)
        zm = lp.zero_mon
        # separate the constant term of the series
        # compute the univariate series f(_x, .., 'x', sum(nv))
        # or _x.f(..., 'x', sum(nv)
        if zm in p:
            x1 = _x + p[zm]
            p1 = p - p[zm]
        else:
            x1 = _x
            p1 = p
        if isinstance(f, str):
            q = getattr(x1, f)(*args1)
        else:
            q = f(x1, *args1)
        a = list(zip(list(q.keys()), list(q.values())))
        a.sort()
        c = [0]*h
        for x in a:
            c[x[0][0]] = x[1]
        p1 = p1.series_from_list(c, args[-2], args[-1])
        return p1

    def _nth_root1(p, n, iv, nv):
        """univariate series nth root of p on commuting ring
        n  integer; compute p^(1/n)
        iv name of the series variable
        nv precision of the series

        The Newton method is used.
        """
        lp = p.lp
        zm = lp.zero_mon
        if zm not in p:
            raise NotImplementedError('no constant term in series')
        if n != int(n):
            raise NotImplementedError('n must be an integer')
        else:
            n = int(n)
        if p[zm] != 1:
            raise NotImplementedError('the constant term must be 1')
        p1 = lp(1)
        if p == 1:
            return p
        if n == 0:
            if p != 0:
                return lp(1)
            else:
                raise ValueError('0^0 expression')
        if n == 1:
            return p
        if n < 0:
            n = -n
            sign = 1
        else:
            sign = 0
        for prec in giant_steps(nv):
            tmp = p1.pow_trunc(n+1, iv, prec)
            tmp = tmp.mul_trunc(p, iv, prec)
            p1 += p1/n - tmp/n
        if sign:
            return p1
        else:
            return p1._series_inversion1(iv, nv)
        return p1

    def sqrt(p, iv, nv):
        """square root of multivariate series p
        iv list of variable names or variable indices
        nv list of truncations for these variables
        In the case of one variable it can also be:
          iv variable name or variable index (0)
          nv truncation integer for the variable
        p is a series with O(x_1^n_1*..x_m^n_m) in
        variables x_k with index or name iv[k-1]
        p has constant term equal to 1

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = (1 + x + x*y).sqrt('x', 3)
        >>> str(p)
        ' -1/8*x^2*y^2 -1/4*x^2*y -1/8*x^2 +1/2*x*y +1/2*x +1'
        """
        p1 = p.nth_root(-2, iv, nv)
        return p.mul_trunc(p1, iv, nv)


    def nth_root(p, n, iv, nv):
        """multivariate series nth root of p
        n  integer; compute p^(1/n)
        iv list of variable names or variable indices
        nv list of truncations for these variables
        In the case of one variable it can also be:
          iv variable name or variable index (0)
          nv truncation integer for the variable
        p is a series with O(x_1^n_1*..x_m^n_m) in
        variables x_k with index or name iv[k-1]
        p has constant term equal to 1
        TODO: case of constant term different from 1

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ, lex)
        >>> p = (1 + x + x*y).nth_root(-3, 'x', 3)
        >>> str(p)
        ' +2/9*x^2*y^2 +4/9*x^2*y +2/9*x^2 -1/3*x*y -1/3*x +1'
        """
        lp = p.lp
        if (p-1).has_constant_term(iv):
            if not lp.SR:
                raise TaylorEvalError('p-1 must not have a constant term in the series variables')
            else:
                if lp.zero_mon in p:
                    c = p[lp.zero_mon]
                    if c.is_positive:
                        return (p/c).nth_root(n, iv, nv)*c**Rational(1, n)
                    else:
                        raise NotImplementedError
        if isinstance(nv,int) and lp.commuting:
            return p._nth_root1(n, iv, nv)
        return p.fun('_nth_root1', n, iv, nv)

    def _log_series0(p, iv, nv):
        p = p - 1
        lp = p.lp
        s = lp(0)
        q = lp(-1)
        n = 1
        z = lp(0)
        while 1:
            q = -q.mul_trunc(p, iv, nv)
            if q == z:
                break
            s += q/n
            n += 1
        return s

    def re_acoth(p, iv, nv):
        """Re(acoth(p)) = 1/2*(log(1+p) - log(1-p))

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = x.re_acoth('x', 8)
        >>> str(p)
        ' +1/7*x^7 +1/5*x^5 +1/3*x^3 +x'
        """
        if p.has_constant_term(iv):
            raise NotImplementedError('p must not have a constant term in the series variables')
        p = p
        lp = p.lp
        s = p.copy()
        q = p
        n = 3
        z = lp(0)
        p2 = p.square_trunc(iv, nv)
        while 1:
            q = q.mul_trunc(p2, iv, nv)
            if q == z:
                break
            s += q/n
            n += 2
        return s

    def acot1(p, iv, nv):
        """acot1(p) = acot(p) - constant part = -p +p^3/3 -p^5/5 + ...
        """
        if p.has_constant_term(iv):
            raise NotImplementedError('p must not have a constant term in the series variables')
        p = p
        lp = p.lp
        s = -p
        q = -p
        n = 3
        z = lp(0)
        p2 = -p.square_trunc(iv, nv)
        while 1:
            q = q.mul_trunc(p2, iv, nv)
            if q == z:
                break
            s += q/n
            n += 2
        return s

    def _log_series(p, iv, nv):
        if len(p) < 50:
            return p._log_series0(iv, nv)
        lp = p.lp
        if lp.base_ring:
            one = lp.base_ring(1)
        else:
            one = lp.ring(1)
        c = [one*0]
        if isinstance(nv,int):
            m = nv
        else:
            m = sum(nv)
        for k in range(1, m):
            if k%2 == 0:
                cf = -one/k
            else:
                cf = one/k
            c.append(cf)
        return (p-1).series_from_list(c, iv, nv)

    def log(p, iv, nv):
        """
        logarithm of p with truncation
        p polynomial starting with 1

        For univariate series or with the total degree
        truncation integral dx p^-1*d p/dx is used.

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = (1+x).log('x', 8)
        >>> str(p)
        ' +1/7*x^7 -1/6*x^6 +1/5*x^5 -1/4*x^4 +1/3*x^3 -1/2*x^2 +x'
        """
        lp = p.lp
        if (p-1).has_constant_term(iv):
            if not lp.SR:
                raise TaylorEvalError('p-1 must not have a constant term in the series variables')
            else:
                if lp.zero_mon in p:
                    c = p[lp.zero_mon]
                    if c.is_positive:
                        return log(c) + (p/c).log(iv, nv)
                    else:
                        raise NotImplementedError
            raise NotImplementedError('p-1 must not have a constant term in the series variables')
        lp = p.lp
        if iv in lp.pol_gens and isinstance(nv,int) and lp.commuting:
            dlog = p.derivative(iv)
            dlog = dlog.mul_trunc(p._series_inversion1(iv, nv), iv, nv-1)
            return dlog.integrate(iv)
        return p._log_series(iv, nv)

    def _atan_series0(p, iv, nv):
        lp = p.lp
        s = lp(1)
        q = lp(1)
        n = 3
        p2 = p.square_trunc(iv, nv)
        while 1:
            q = -q.mul_trunc(p2, iv, nv)
            if not q:
                break
            s += q/n
            n += 2
        s = s.mul_trunc(p, iv, nv)
        return s

    def _atan_series(p, iv, nv):
        if len(p) < 50:
            return p._atan_series0(iv, nv)
        lp = p.lp
        mo = lp.ring(-1)
        c = [-mo]
        if isinstance(nv,int):
            m = nv
        else:
            m = sum(nv)
        p2 = p.square_trunc(iv, nv)
        for k in range(1, m):
            c.append(mo**k/(2*k+1))
        s = p2.series_from_list(c, iv, nv)
        s = s.mul_trunc(p, iv, nv)
        #s += p.trunc(iv, nv)
        return s


    def atan(p, iv, nv):
        """
        arctangent of a series
        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = x.atan('x', 8)
        >>> str(p)
        ' -1/7*x^7 +1/5*x^5 -1/3*x^3 +x'
        """
        if p.has_constant_term(iv):
            raise NotImplementedError('polynomial must not have constant term in the series variables')
        lp = p.lp
        if iv in lp.pol_gens and isinstance(nv,int) and lp.commuting:
            dp = p.derivative(iv)
            p1 = p.square_trunc(iv, nv) + 1
            p1 = p1.series_inversion(iv, nv-1)
            p1 = dp.mul_trunc(p1, iv, nv-1)
            return p1.integrate(iv)
        else:
            return p._atan_series(iv, nv)

    def lambert(p, iv, nv):
        """
        principal branch of the Lambert W function

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = x.lambert('x', 8)
        >>> str(p)
        ' +16807/720*x^7 -54/5*x^6 +125/24*x^5 -8/3*x^4 +3/2*x^3 -x^2 +x'
        """
        lp = p.lp
        p1 = lp(0)
        if p.has_constant_term(iv):
            raise NotImplementedError('polynomial must not have constant term in the series variables')
        if iv in lp.pol_gens and isinstance(nv,int) and lp.commuting:
            for prec in giant_steps(nv):
                e = p1.exp(iv, prec)
                p2 = e.mul_trunc(p1, iv, prec) - p
                p3 = e.mul_trunc(p1 + 1, iv, prec)
                p3 = p3.series_inversion(iv, prec)
                tmp = p2.mul_trunc(p3, iv, prec)
                p1 -= tmp
            return p1
        else:
            raise NotImplementedError

    def asin(p, iv, nv):
        """
        arcsin of a series

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = x.asin('x', 8)
        >>> str(p)
        ' +5/112*x^7 +3/40*x^5 +1/6*x^3 +x'
        """
        if p.has_constant_term(iv):
            raise NotImplementedError('polynomial must not have constant term in the series variables')
        lp = p.lp
        if iv in lp.pol_gens and isinstance(nv,int) and lp.commuting:
            # get a good value
            if len(p) > 20:
                dp = p.derivative(iv)
                p1 = 1 - p.square_trunc(iv, nv-1)
                p1 = p1.nth_root(-2, iv, nv-1)
                p1 = dp.mul_trunc(p1, iv, nv-1)
                return p1.integrate(iv)
            if lp.base_ring:
                one = lp.base_ring(1)
            else:
                one = lp.ring(1)
            c = [0, one, 0]
            if isinstance(nv,int):
                m = nv
            else:
                m = sum(nv)

            for k in range(3, m, 2):
                if k%2 == 1:
                    c.append((k-2)**2*c[-2]/(k*(k-1)))
                    c.append(0)
            return p.series_from_list(c, iv, nv)

        else:
            raise NotImplementedError

    def asinh(p, iv, nv):
        """
        arcsinh of a series

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = x.asinh('x', 8)
        >>> str(p)
        ' -5/112*x^7 +3/40*x^5 -1/6*x^3 +x'
        """
        if p.has_constant_term(iv):
            raise NotImplementedError('polynomial must not have constant term in the series variables')
        lp = p.lp
        if iv in lp.pol_gens and isinstance(nv,int) and lp.commuting:
            dp = p.derivative(iv)
            p1 = 1 + p.square_trunc(iv, nv-1)
            p1 = p1.nth_root(-2, iv, nv-1)
            p1 = dp.mul_trunc(p1, iv, nv-1)
            return p1.integrate(iv)
        else:
            raise NotImplementedError

    def _atanh_series0(p, iv, nv):
        lp = p.lp
        s = lp(1)
        q = lp(1)
        n = 3
        p2 = p.square_trunc(iv, nv)
        while 1:
            q = q.mul_trunc(p2, iv, nv)
            if not q:
                break
            s += q/n
            n += 2
        s = s.mul_trunc(p, iv, nv)
        return s

    def _atanh_series(p, iv, nv):
        if len(p) < 50:
            return p._atanh_series0(iv, nv)
        lp = p.lp
        one = lp.ring(1)
        c = [one]
        if isinstance(nv,int):
            m = nv
        else:
            m = sum(nv)
        p2 = p.square_trunc(iv, nv)
        for k in range(1, m):
            c.append(one/(2*k+1))
        s = p2.series_from_list(c, iv, nv)
        s = s.mul_trunc(p, iv, nv)
        return s

    def atanh(p, iv, nv):
        """ hyperbolic arctangent

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = x.atanh('x', 8)
        >>> str(p)
        ' +1/7*x^7 +1/5*x^5 +1/3*x^3 +x'
        """
        if p.has_constant_term(iv):
            raise NotImplementedError('polynomial must not have constant term in the series variables')
        lp = p.lp
        if iv in lp.pol_gens and isinstance(nv,int) and lp.commuting:
            dp = p.derivative(iv)
            p1 = -p.square_trunc(iv, nv) + 1
            p1 = p1.series_inversion(iv, nv-1)
            p1 = dp.mul_trunc(p1, iv, nv-1)
            return p1.integrate(iv)
        else:
            return p._atanh_series(iv, nv)

    def _tanh1(p, iv, nv):
        lp = p.lp
        p1 = lp(0)
        for prec in giant_steps(nv):
            tmp = p - p1.atanh(iv, prec)
            tmp = tmp.mul_trunc(1 - p1.square(), iv, prec)
            p1 += tmp
        return p1

    def tanh(p, iv, nv):
        """ hyperbolic tangent of a series

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = x.tanh('x', 8)
        >>> str(p)
        ' -17/315*x^7 +2/15*x^5 -1/3*x^3 +x'
        """
        lp = p.lp
        if p.has_constant_term(iv):
            raise NotImplementedError('p must not have constant part in series variables')
        p1 = lp(0)
        if isinstance(nv,int) and lp.commuting:
            if iv in lp.pol_gens:
                return p._tanh1(iv, nv)
        if isinstance(nv,int):
            nv = [int(nv)]
        return p.fun('tanh', iv, nv)

    def _tan1(p, iv, nv):
        lp = p.lp
        p1 = lp(0)
        for prec in giant_steps(nv):
            tmp = p - p1.atan(iv, prec)
            tmp = tmp.mul_trunc(1 + p1.square(), iv, prec)
            p1 += tmp
        return p1

    def tan(p, iv, nv):
        """tangent of a series

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = x.tan('x', 8)
        >>> str(p)
        ' +17/315*x^7 +2/15*x^5 +1/3*x^3 +x'
        """
        lp = p.lp
        if p.has_constant_term(iv):
            raise NotImplementedError('p must not have constant part in series variables')
        p1 = lp(0)
        if isinstance(nv,int) and lp.commuting:
            if iv in lp.pol_gens and isinstance(nv,int):
                return p._tan1(iv, nv)
        if isinstance(nv,int):
            nv = [int(nv)]
        return p.fun('tan', iv, nv)

    def _exp_series0(p, iv, nv):
        if 0 in p:
            raise NotImplementedError('p must not have constant part')
        lp = p.lp
        q = lp(1)
        n = 1
        k = 1
        p1 = lp(1)
        z = lp(0)
        while 1:
            q = q.mul_trunc(p, iv, nv)
            if not q:
                break
            p1 += q/n
            k += 1
            n *= k
        return p1


    def _exp_series(p, iv, nv):
        lp = p.lp
        zm = lp.zero_mon
        if zm in p:
            raise NotImplementedError('p must not have constant part')
        n = lp.ring(1)
        k = 1
        c = []
        if isinstance(nv,int):
            m = nv
        else:
            m = sum(nv)
        for k in range(m):
            c.append(1/n)
            k += 1
            n *= k
        return p.series_from_list(c, iv, nv)


    def _exp1(p, iv, nv):
        """
        exponential of a univariate series, or of a multivariate
        series with total degree truncation
        """
        lp = p.lp
        zm = lp.zero_mon
        # TODO p must not have terms not containing series variables
        if zm in p:
            raise NotImplementedError('p must not have constant part')
        p1 = lp(1)
        for prec in giant_steps(nv):
            tmp = (p - p1.log(iv, prec)).mul_trunc(p1, iv, prec)
            p1 += tmp
        return p1

    def exp(p, iv, nv):
        """
        exponential of a series

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = x.exp('x', 6)
        >>> str(p)
        ' +1/120*x^5 +1/24*x^4 +1/6*x^3 +1/2*x^2 +x +1'
        """
        lp = p.lp
        if p.has_constant_term(iv):
            zm = lp.zero_mon
            if not lp.SR:
                raise TaylorEvalError('p must not have constant part in series variables')
            return exp(p[zm])*(p - p[zm]).exp(iv, nv)
        p1 = lp(1)
        if isinstance(nv,int) and len(p) > 20 and lp.commuting:
            if iv in lp.pol_gens and isinstance(nv,int):
                return p._exp1(iv, nv)
        if lp.base_ring:
            one = lp.base_ring(1)
        else:
            one = lp.ring(1)
        n = 1
        k = 1
        c = []
        if isinstance(nv,int):
            m = nv
        else:
            m = sum(nv)
        for k in range(m):
            c.append(one/n)
            k += 1
            n *= k
        return p.series_from_list(c, iv, nv)

    def sin(p, iv, nv):
        """ sin of a series

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = x.sin('x', 6)
        >>> str(p)
        ' +1/120*x^5 -1/6*x^3 +x'
        """
        lp = p.lp
        if p.has_constant_term(iv):
            zm = lp.zero_mon
            if not lp.SR:
                raise TaylorEvalError
            c = p[zm]
            if c.is_number and not c.is_real:
                raise TaylorEvalError
            p1 = p - c
            return cos(c)*p1.sin(iv, nv) + sin(c)*p1.cos(iv, nv)
        # get a good value
        if len(p) > 20:
            t = (p/2).tan(iv, nv)
            t2 = t.square_trunc(iv, nv)
            p1 = (1+t2).series_inversion(iv, nv)
            return p1.mul_trunc(2*t, iv, nv)
        if lp.base_ring:
            one = lp.base_ring(1)
        else:
            one = lp.ring(1)
        n = 1
        c = [0]
        if isinstance(nv,int):
            m = nv
        else:
            m = sum(nv)
        for k in range(2, m+2, 2):
            c.append(one/n)
            c.append(0)
            n *= -k*(k+1)
        return p.series_from_list(c, iv, nv)

    def cos(p, iv, nv):
        """ cos of a series

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = x.cos('x', 6)
        >>> str(p)
        ' +1/24*x^4 -1/2*x^2 +1'
        """
        lp = p.lp
        if p.has_constant_term(iv):
            zm = lp.zero_mon
            if not lp.SR:
                raise TaylorEvalError
            c = p[zm]
            if not c.is_real:
                raise NotImplementedError
            p1 = p - c
            return cos(c)*p1.cos(iv, nv) - sin(c)*p1.sin(iv, nv)
        # get a good value
        if len(p) > 20:
            t = (p/2).tan(iv, nv)
            t2 = t.square_trunc(iv, nv)
            p1 = (1+t2).series_inversion(iv, nv)
            return p1.mul_trunc(1-t2, iv, nv)
        if lp.base_ring:
            one = lp.base_ring(1)
        else:
            one = lp.ring(1)
        n = 1
        c = []
        if isinstance(nv,int):
            m = nv
        else:
            m = sum(nv)
        for k in range(2, m+2, 2):
            c.append(one/n)
            c.append(0)
            n *= -k*(k-1)
        return p.series_from_list(c, iv, nv)


    def cos_sin(p, iv, nv):
        """
        tuple (p.cos(iv, iv), p.sin(iv, iv))
        """
        t = (p/2).tan(iv, nv)
        t2 = t.square_trunc(iv, nv)
        p1 = (1+t2).series_inversion(iv, nv)
        return (p1.mul_trunc(1-t2, iv, nv), p1.mul_trunc(2*t, iv, nv))

    def sinh(p, iv, nv):
        """ hyperbolic sin of a series

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = x.sinh('x', 8)
        >>> str(p)
        ' +1/5040*x^7 +1/120*x^5 +1/6*x^3 +x'
        """
        t = p.exp(iv, nv)
        t1 = t.series_inversion(iv, nv)
        return (t - t1)/2

    def cosh(p, iv, nv):
        """ cos of a series

        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ, lex)
        >>> p = x.cosh('x', 8)
        >>> str(p)
        ' +1/720*x^6 +1/24*x^4 +1/2*x^2 +1'
        """
        t = p.exp(iv, nv)
        t1 = t.series_inversion(iv, nv)
        return (t + t1)/2

    def cosh_sinh(p, iv, nv):
        """ tuple (p.cosh(iv, iv), p.sinh(iv, iv))
        """
        t = p.exp(iv, nv)
        t1 = t.series_inversion(iv, nv)
        return (t + t1)/2, (t - t1)/2

    def toSR(p, lp1):
        """convert coefficients to Rational
        """
        assert lp1.SR
        lp = p.lp
        if lp == lp1:
            return p
        p1 = lp1(0)
        for expv, c in p.items():
            p1[expv] = QQ.to_sympy(c)
        return p1


    def tobasic(p, *gens):
        """Convert a Poly into a Sympy expression.
        >>> from sympy.core.numbers import Rational as QQ
        >>> from sympy.polys.monomialtools import lex
        >>> from sympy.polys.lpoly import lgens
        >>> from sympy import Symbol
        >>> lp, X = lgens('X', QQ, lex)
        >>> x = Symbol('x')
        >>> p = (1 + X)**3
        >>> str(p)
        ' +X^3 +3*X^2 +3*X +1'
        >>> p1 = p.tobasic(x)
        >>> p1
        x**3 + 3*x**2 + 3*x + 1
        """
        ring = p.lp.ring
        if str(ring) == 'QQ':
            def convert(c):
                return ring.to_sympy(c)
        else:
            def convert(c):
                return c
        result = []
        for monom, coeff in p.items():
            term = [convert(coeff)]
            for g, m in zip(gens, monom):
                term.append(Pow(g, m))

            result.append(Mul(*term))
        return Add(*result)

class LPolySubs(object):
    """class for substitutions of variables with polynomials,
    possibly truncated.

    >>> from sympy.core.numbers import Rational as QQ
    >>> from sympy.polys.monomialtools import lex
    >>> from sympy.polys.lpoly import lgens
    >>> from sympy.polys.lpoly import LPolySubs
    >>> lp, x, y = lgens('x, y', QQ, lex)
    >>> sb = LPolySubs(lp, lp, {'y':y+1})
    >>> p1 = sb.subs((x + y)**2)
    >>> str(p1)
    ' +x^2 +2*x*y +2*x +y^2 +2*y +1'
    >>> p1 = sb.subs_trunc((x + y)**4, 'x', 2)
    >>> str(p1)
    ' +4*x*y^3 +12*x*y^2 +12*x*y +4*x +y^4 +4*y^3 +6*y^2 +4*y +1'
    """
    def __init__(self, lp1, lp2, rules):
        self.lp1 = lp1
        gens_dict = lp1.gens_dict
        self.lp2 = lp2
        if lp1.ring != lp2.ring:
            raise NotImplementedError
        d = {} # replace monomials with (i, pw)
        gens = lp1.gens()
        for i in range(lp1.ngens):
            d[(i, 1)] = gens[i]
        for var in rules:
            d[(gens_dict[var], 1)] = rules[var]
        self.d = d

    def subs(self, p):
        lp1 = self.lp1
        lp2 = self.lp2
        ngens = lp1.ngens
        assert p.lp == lp1
        d = self.d.copy()
        p1 = Poly(lp2)
        for expv in p:
            p2 = lp2(1)
            for i in range(ngens):
                pw = expv[i]
                if pw == 0:
                    continue
                if (i, pw) not in d:
                    if pw%2 == 0 and (i, pw/2) in d:
                        d[(i, pw)] = d[(i, pw/2)]**2
                    elif (i, pw-1) in d:
                        d[(i, pw)] = d[(i, pw-1)]*d[(i, 1)]
                    else:
                        d[(i, pw)] = d[(i, 1)]**pw
                p2 *= d[(i, pw)]
            p1 += p2*p[expv]
        return p1

    def subs_trunc(self, p, ii, h):
        """substitution with truncation of variable(s) corresponding
        to ii and truncation order(s) h
        """
        lp1 = self.lp1
        lp2 = self.lp2
        ngens = lp1.ngens
        assert p.lp == self.lp1
        d = self.d.copy()
        p1 = lp2(0)
        for expv in p:
            p2 = lp2(1)
            for i in range(ngens):
                pw = expv[i]
                if pw == 0:
                    continue
                if (i, pw) not in d:
                    if pw%2 == 0 and (i, pw/2) in d:
                        d[(i, pw)] = d[(i, pw/2)].square_trunc(ii, h)
                    elif (i, pw-1) in d:
                        d[(i, pw)] = d[(i, pw-1)].mul_trunc(d[(i, 1)], ii, h)
                    else:
                        d[(i, pw)] = d[(i, 1)].pow_trunc(pw, ii, h)
                p2 = p2.mul_trunc(d[(i, pw)], ii, h)
            p1 += p2*p[expv]
        return p1
