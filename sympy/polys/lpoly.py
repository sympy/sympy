""" lpoly translation of rmpoly, see http://code.google.com/p/rmpoly,
using monomials instead of packed exponents for the monomials.
"""

import sympy
from sympy import S, Rational, sign, symbols, Add, Mul, Pow, Symbol
from sympy.core.compatibility import is_sequence
from sympy.polys.monomialtools import monomial_mul, monomial_div, monomial_min, lex
from sympy.polys.domains import ZZ, QQ, PythonRationalType
from sympy.ntheory.residue_ntheory import int_tested

from copy import copy
import re
import math

_rpm = re.compile('[+-]')
_re_var_split = re.compile(r"\s*,\s*|\s+")

class TaylorEvalError(TypeError):
    """
    Exception used in ltaylor.taylor
    """
    pass

def giant_steps(target):
    """
    list of precision steps for the Newton's method

    code adapted from mpmath/libmp/libintmath.py
    """
    L = [target]
    start = 2
    while 1:
        Li = L[-1]//2 + 2
        if Li >= L[-1] or Li < start:
            if L[-1] != start:
                L.append(start)
            break
        L.append(Li)
    return L[::-1]

def monomial_zero(n):
    """zero monomial in n variables"""
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

def monomial_as_expr(monom, *gens):
    """
    >>> from sympy.polys.lpoly import monomial_as_expr
    >>> from sympy import symbols
    >>> x, y = symbols('x, y')
    >>> monomial_as_expr((2, 1), x, y)
    x**2*y
    """
    assert len(monom) == len(gens)
    return Mul(*[Pow(g, m) for (g, m) in zip(gens, monom)])

gmpy_mode = not QQ(1).__class__ is PythonRationalType

def PythonRationalType_new(p):
    """
    returns a PythonRationalType
    if p is a PythonRationalType it returns it
    if p is a sring representing a rational number or an integer,
    it returns a PythonRationalType object

    Examples
    ========
    >>> from sympy.polys.domains.pythonrationalfield import PythonRationalType
    >>> from sympy.polys.lpoly import PythonRationalType_new
    >>> a =  PythonRationalType(2, 3)
    >>> PythonRationalType_new(a)
    2/3
    """
    if isinstance(p, PythonRationalType):
        return p
    elif isinstance(p, basestring):
        a = p.split('/')
        if len(a) == 2:
            return PythonRationalType(int(a[0]), int(a[1]))
        else:
            return PythonRationalType(int(a[0]), 1)
    else:
        return PythonRationalType(p)

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
    construct elements of LPolyElement, a multivariate polynomial ring with
    monomial ordering.
    lp.ring generates the elements in K

      lp.ring(0)  zero element in K
      lp.ring(1)  unit element in K

    O order object
    The string representation of polynomials is in decreasing order.

    Examples
    ========
    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import LPoly
    >>> lp = LPoly('x, y', QQ)
    >>> lp.zero_mon
    (0, 0)
    >>> x, y = lp.gens
    >>> p = (x + y)**2
    >>> p
    x**2 + 2*x*y + y**2
    """

    def __init__(self, pol_gens, ring, order, **kwds):
        if not is_sequence(pol_gens, include=(str, Symbol)) or not pol_gens:
            raise ValueError('expecting string, Symbol or other ordered iterable')

        if isinstance(pol_gens, str):
            pol_gens = [s.name for s in symbols(pol_gens, seq=True)]
        elif isinstance(pol_gens, Symbol):
            pol_gens = pol_gens.name
        elif isinstance(pol_gens[0], Symbol):
            pol_gens = [s.name for s in pol_gens]
        elif not isinstance(pol_gens[0], str):
            raise ValueError('expecting iterables of Symbols or strings')

        self.pol_gens = tuple(pol_gens)
        self.ngens = len(pol_gens)
        self.ring = ring
        self.order = order
        self.gens_dict = dict(zip(self.pol_gens, xrange(self.ngens)))
        self.parens = kwds.get('parens', False)
        self.commuting = kwds.get('commuting', True)
        if isinstance(ring, BaseLPoly):
            self.parens = True
            if not ring.commuting:
                self.commuting = False
        str_ring = str(self.ring)
        if 'mpc' in str_ring or 'complex' in str_ring.lower():
            self.parens = True
        self.SR = False
        if hasattr(ring, '__name__'):
            if ring.__name__ == 'sympify':
                self.SR = True
                self.parens = True
        self.ring_new = ring
        if not self.SR and not gmpy_mode:
            if ring.__class__ == QQ.__class__ and ring == QQ:
                self.ring_new = PythonRationalType_new
        self.zero_mon = monomial_zero(self.ngens)

    def __str__(self):
        try:
            ring_name = self.ring.__name__
        except AttributeError:
            ring_name = str(self.ring)
        s = 'LPoly with ngens=%d ring=%s' % (self.ngens, ring_name)

        return s

    @property
    def gens(self):
        """return a list of the LPoly generators
        """
        one = self.ring(1)
        ngens = self.ngens
        a = []
        for i in range(ngens):
            expv = monomial_basis(i, ngens)
            p = LPolyElement(self)
            p[expv] = one
            a.append(p)
        return tuple(a)

    def read_monom(self, s):
        """compute expv for the string 'x_0**e_0*...'
        where expv is (e_0, ..)
        """
        gens_dict = self.gens_dict
        ngens = self.ngens
        s = s.replace(' ', '')
        s = s.replace("**", "^")
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
        """compute the tuple (expv, c) for a string 'c*x_0**e_0*...'
        where expv = (e0, ..)
        """
        gens_dict = self.gens_dict
        ngens = self.ngens
        ring_new = self.ring_new
        s = s.replace(' ', '')
        s = s.replace('**', '^')
        if s[0].isdigit() and '*' not in s:
            return (self.zero_mon, ring_new(s))
        if s[0].isdigit():
            a = s.split('*', 1)
            coeff = ring_new(a[0])
            a = a[1].split('*')
        else:
            a = s.split('*')
            coeff = self.ring(1)
        expv = [0]*ngens
        # e.g. a = ['x0**6', 'x1**6', 'x10**2']
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

    def __call__(self, s):
        """
        Generate a polynomial from a string.

        The string must be in the form
            c_0*m_0 + c_1*m_1 + ...
        where c_i are coefficients and m_i are monomials

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = lp('1/2 + 2/3*x**4 + 1/5*x*y**2')
        >>> p
        2/3*x**4 + 1/5*x*y**2 + 1/2
        """
        if isinstance(s, LPolyElement):
            if s.lp == self:
                return s
            else:
                raise NotImplementedError('cannot convert polynomial')
        if not isinstance(s, str):
            c = self.ring_new(s)
            p = LPolyElement(self)
            if c:
                p[self.zero_mon] = c
            return p
        if self.parens:
            raise NotImplementedError
        s = s.replace(' ', '')
        a = []
        if s[0] in '+-':
            sgn = s[0]
        else:
            sgn = '+'
            s = sgn + s
        it = _rpm.finditer(s)
        next(it)
        s1 = s[1:]
        if _rpm.search(s1):
            start = 1
            for match in it:
                t = match.span()
                a.append((sgn, s[start:t[0]]))
                start = t[0] + 1
                sgn = s[t[0]]
            a.append((sgn, s[t[1]:]))
        else:
            a.append((sgn, s1))
        p = LPolyElement(self)
        for sgn, s1 in a:
            m = self.mon_eval(s1)
            if sgn == '-':
                m = (m[0], -m[1])
            if m[0] in p:
                p[m[0]] += m[1]
            else:
                p[m[0]] = m[1]
        p.strip_zero()
        return p

    def from_mon(self, a):
        """polynomial from the monomial coeff*x_i**j
           a = (i, j, coeff)
        """
        p = LPolyElement(self)
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

    def from_dict(self, d):
        """
        polynomial from a dictionary

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> from sympy.abc import x, y
        >>> p = (x + x*y + x**2*y + x**3*y**2).as_poly()
        >>> d = p.as_dict()
        >>> lp, X, Y = lgens('X, Y', QQ)
        >>> lp.from_dict(d)
        X**3*Y**2 + X**2*Y + X*Y + X
        """
        p = LPolyElement(self)
        for k, v in d.iteritems():
            if self.SR:
                p[k] = v
            elif isinstance(v, Rational):
                v = QQ(v.p, v.q)
            else:
                raise TaylorEvalError
            p[k] = v
        p.strip_zero()
        return p

class LPoly(BaseLPoly):
    """class of polynomials on a ring with elements with one term

    Examples
    ========
    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens
    >>> lp, x, y = lgens('x, y', QQ)
    >>> lp.__class__
    <class 'sympy.polys.lpoly.LPoly'>
    """
    def __init__(self, pol_gens, ring, order=lex, **kwds):
        BaseLPoly.__init__(self, pol_gens, ring, order, **kwds)
        self.__name__ = 'LPoly'


class MLPoly(BaseLPoly):
    """class of polynomials on a ring with elements with more than one term

    This class differs from LPoly only because the coefficients of
    its polynomials are surrounded by parenthesis

    Examples
    ========
    >>> from sympy.polys.lpoly import lgens, MLPoly
    >>> from sympy.core import sympify
    >>> from sympy.functions.elementary.exponential import log
    >>> lp, x, y = lgens('x, y', sympify)
    >>> p = (x + x*log(2) + y)**2
    >>> lp.__class__
     <class 'sympy.polys.lpoly.MLPoly'>
    """
    def __init__(self, pol_gens, ring, order, **kwds):
        BaseLPoly.__init__(self, pol_gens, ring, order, **kwds)
        self.parens = True
        self.__name__ = 'MLPoly'


class NCLPoly(BaseLPoly):
    """class of polynomials on noncommuting ring
    """
    def __init__(self, pol_gens, ring, order, **kwds):
        BaseLPoly.__init__(self, pol_gens, ring, order, **kwds)
        self.parens = True
        self.commuting = False
        self.__name__ = 'NCLPoly'

def lgens(pol_gens, ring, order=lex, **kwds):
    """
     factory function to generate LPoly object and its generators

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens
    >>> lp, x, y = lgens('x, y', QQ)
    >>> (x + y)**2
    x**2 + 2*x*y + y**2
    """
    lp = BaseLPoly(pol_gens, ring, order, **kwds)
    if lp.parens:
        lp.__class__ = MLPoly
        lp.__name__ = 'MLPoly'
    else:
        lp.__class__ = LPoly
        lp.__name__ = 'LPoly'
    return (lp,) + lp.gens

def mlgens(pol_gens, ring, order=lex, **kwds):
    """
    factory function to generate MLPoly object and its generators
    """
    lp = MLPoly(pol_gens, ring, order, **kwds)
    return (lp,) + lp.gens

def nclgens(pol_gens, ring, order=lex, **kwds):
    """
    factory function to generate NCLPoly object and its generators
    """
    lp = NCLPoly(pol_gens, ring, order, **kwds)
    return (lp,) + lp.gens




class LPolyElement(dict):
    """
    elements of a multivariate polynomial ring

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens, LPolyElement
    >>> lp, x, y = lgens('x, y', QQ)
    >>> p = LPolyElement(lp)
    >>> p
    0
    >>> p1 = x + y**2
    >>> isinstance(p1, LPolyElement)
    True
    """
    def __init__(self, lp, **kw):
        self.lp = lp
        dict.__init__(self, **kw)

    def __str__(self):
        """string representation of a polynomial
        Terms are in decreasing order; term coefficients are in front

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x**2/3 + 2*x*y**2/5
        >>> p
        1/3*x**2 + 2/5*x*y**2
        """
        if not self:
            return '0'
        lp = self.lp
        if lp.parens:
            return self._tostr()
        pol_gens = lp.pol_gens
        ngens = lp.ngens
        zm = self.lp.zero_mon
        s = []
        a = list(self.keys())
        a.sort(key=lambda m: lp.order(m), reverse=True)
        for expv in a:
            c = self[expv]
            if c > 0:
                s += ' + '
            else:
                s += ' - '
            if c < 0:
                c = -c
            if c != 1 or expv == zm:
                cnt1 = self.__coefficient_to_str(c)
            else:
                cnt1 = ''
            sa = []
            for i in range(ngens):
                exp = expv[i]
                if exp > 1:
                    sa.append('%s**%d' % (pol_gens[i], exp))
                if exp == 1:
                    sa.append('%s' % pol_gens[i])
            if cnt1:
                sa = [cnt1] + sa
            s += '*'.join(sa)
        # remove leading ' + ', and correct ' - ' to '-'
        s = ''.join(s)
        if s[:3] == " + ":
            s = s[3:]
        elif s[:3] == " - ":
            s = "-" + s[3:]
        return s

    __repr__ = __str__

    def __coefficient_to_str(self, c):
        """
        Represent coefficient as string.
        """
        if isinstance(c, PythonRationalType) and c.q == 1:
            return str(c.p)

        return str(c)

    def _tostr(self):
        lp = self.lp
        pol_gens = lp.pol_gens
        ngens = lp.ngens
        zm = lp.zero_mon
        s = []
        a = list(self.keys())
        a.sort(key=lambda m: lp.order(m), reverse=True)
        for expv in a:
            c = self[expv]
            s.append(' + (%s)' % c)
            if expv != zm:
                s += '*'
            # throw away the bits for the hidden variable
            i = 0
            sa = []
            for i in range(ngens):
                exp = expv[i]
                if exp > 1:
                    sa.append('%s**%d' % (pol_gens[i], exp))
                if exp == 1:
                    sa.append('%s' % pol_gens[i])
            s.append('*'.join(sa))
        s = ''.join(s)
        s = s.lstrip(' +')
        return s

    def strip_zero(p):
        """eliminate monomials with zero coefficient
        """
        z = p.lp.ring(0)
        for k, v in list(p.items()):
            if v == z:
                del p[k]

    def variables(self):
        """return the tuple of the indices of the
        variables occurring in the polynomial p

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y, z = lgens('x, y, z', QQ)
        >>> p = x + y**2 + x*y
        >>> p.variables()
        (0, 1)
        """
        a = []
        for expv in self:
            i = 0
            for i, e in enumerate(expv):
                if e and i not in a:
                    a.append(i)
        a.sort()
        return tuple(a)

    def __eq__(p1, p2):
        """
        equality test for polynomials

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p1 = (x + y)**2 + (x - y)**2
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
        if isinstance(p2, LPolyElement):
            if lp1 == p2.lp:
                return dict.__eq__(p1, p2) #or dict.__eq__(p2, p1)
        else:
            zm = lp1.zero_mon
            # assume p2 is a coefficient
            if zm not in p1 or len(p1) > 1:
                return False
            return p1[zm] == lp1.ring_new(p2)

    def __add__(p1, p2):
        """add two polynomials

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = (x + y)**2 + (x - y)**2
        >>> p
        2*x**2 + 2*y**2
        """
        if not p2:
            return p1.copy()
        lp1 = p1.lp
        zm = lp1.zero_mon
        if isinstance(p2, LPolyElement):
            if lp1 == p2.lp:
                p = LPolyElement(lp1)
                for k, v in p1.iteritems():
                    if k in p2:
                        r = v + p2[k]
                        if r:
                            p[k] = r
                    else:
                        p[k] = v
                for k, v in p2.iteritems():
                    if k not in p1:
                        p[k] = v
                return p
            elif p2.lp.__class__ == lp1.ring.__class__ and p2.lp == lp1.ring:
                p = p1.copy()
                if zm not in list(p1.keys()):
                    p[zm] = lp1.ring_new(p2)
                else:
                    if p2 == -p[zm]:
                        del p[zm]
                    else:
                        p[zm] += p2
                return p
            elif lp1.__class__ == p2.lp.ring.__class__ and lp1 == p2.lp.ring:
                return p2 + p1
            else:
                raise ValueError('cannot sum p1 and p2')
        # assume p2 in a number
        else:
            p = p1.copy()
            cp2 = lp1.ring_new(p2)
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
        return a copy of polynomial self

        polynomials are mutable; if one is interested in preserving
        a polynomial, and one plans to use inplace operations, one
        can copy the polynomial

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = (x + y)**2
        >>> p1 = p.copy()
        >>> p2 = p
        >>> p[lp.zero_mon] = 3
        >>> p
        x**2 + 2*x*y + y**2 + 3
        >>> p1
        x**2 + 2*x*y + y**2
        >>> p2
        x**2 + 2*x*y + y**2 + 3
        """
        return copy(self)

    def coefficient(self, p1):
        """the coefficient of a monomial p1 is the sum of the terms in
        self which have the same degrees as the variables present in p1,
        divided by p1

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = (1 + x + y)**3
        >>> p1 = p.coefficient(y**2)
        >>> p1
        3*x + 3
        """
        lp = self.lp
        k = list(p1.keys())
        if len(k) != 1:
            raise TypeError('the argument of coeff must be a monomial')
        expv1 = k[0]
        v1 = p1.variables()
        p = LPolyElement(lp)
        for expv in self:
            for i in v1:
                if expv[i] != expv1[i]:
                    break
            else:
                p[monomial_div(expv, expv1)] = self[expv]
        return p

    def coefficient_t(self, t):
        """coefficient x_i^j of p
        for j=0 it gives the terms independent from the variable x_i

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = (1 + x + y)**3
        >>> p.coefficient_t((0, 0))
        y**3 + 3*y**2 + 3*y + 1
        >>> p.coefficient_t((0, 2))
        3*y + 3
        """
        i, j = t
        expv1 = [0]*self.lp.ngens
        expv1[i] = j
        expv1 = monomial_from_sequence(expv1)
        lp = self.lp
        p = LPolyElement(lp)
        for expv in self:
            if expv[i] == j:
                p[monomial_div(expv, expv1)] = self[expv]
        return p


    def coeff(self, p1):
        """monomial coefficient in self of the leading monomial of p1

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = (1 + x + y)**3/5
        >>> p.coeff(y**2)
        3/5
        >>> p.coeff(x*y)
        6/5
        """
        return self.get(p1.leading_expv(), self.lp.ring(0))

    def series_reversion(p, i, n, j):
        """reversion of a series

        p is a series with O(x**n) of the form p = a*x + f(x)
        where `a` is a number different from 0

        f(x) = sum( a_k*x_k, k in range(2, n))

          a_k can depend polynomially on other variables, not indicated.
          x variable with index i, or with name i
          y variable with index j, or with name j

        solve p = y, that is given
        a*x + f(x) - y = 0
        find the solution x = r(y) up to O(y**n)

        Algorithm:
        if r_i is the solution at order i
        a*r_i + f(r_i) - y = O(y^(i + 1))
        r_(i + 1) = r_i + e such that
        a*r_(i + 1) + f(r_(i + 1)) - y = O(y^(i + 2))
        a*e + f(r_i) = O(y^(i + 2))
        e = -f(r_i)/a
        so that one uses the recursion relation
        r_(i + 1) = r_i -f(r_i)/a
        with the boundary condition
        r_1 = y

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x + x**2
        >>> p1 = p.series_reversion('x', 4, 'y')
        >>> p1
        2*y**3 - y**2 + y
        >>> (p1 + p1**2).trunc('y', 4)
        y
        """
        lp = p.lp
        if not lp.commuting:
            raise NotImplementedError
        if i.__class__ == str:
            nx = i
            i = p.lp.gens_dict[i]
        else:
            nx = lp.pol_gens[i]
        y = lp(j)
        if j.__class__ == str:
            ny = j
            j = p.lp.gens_dict[j]
        else:
            ny = lp.pol_gens[j]
        if p.coefficient_t((i, 0)):
            raise ValueError('part independent from %s must be 0' % nx)
        a = p.coefficient_t((i, 1))
        zm = lp.zero_mon
        assert zm in a and len(a) == 1
        a = a[zm]
        r = lp(ny)/a
        for i in range(2, n):
            sb = LPolySubs(lp, lp, {nx:r})
            sp = sb.subs_trunc(p, ny, i + 1)
            sp = sp.coefficient_t((j, i))*y**i
            r -= sp/a
        return r

    def subs(p, **rules):
        """
        substitution

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x**2 + y**2
        >>> p1 = p.subs(x=x + y, y=x - y)
        >>> p1
        2*x**2 + 2*y**2
        """
        lp = p.lp
        sb = LPolySubs(lp, lp, rules)
        return sb.subs(p)

    def subs_trunc(p, iv, nv, **rules):
        """
        substitution with truncation

          p     input polynomial
          iv    (name of the) variable in which the series truncation is done
          nv    order of the truncation
          rules rules of substitution x=px, y=py, .. where x, y are
                in p.lp.gens and px, py are the polynomials with which
                they are substituted:
                p = ... + c_{i,j}*x**i*y**j + ... ->
                    ... + c_{i,j}*px**i*py**j + ...

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x**2 + y**2
        >>> p2 = p.subs_trunc('x', 3, y=(x + y)**2)
        >>> p2
        6*x**2*y**2 + x**2 + 4*x*y**3 + y**4
        >>> p.subs_trunc('x', 3, y=p.lp(1))
        x**2 + 1

        One can even substitute the series variable with a constant;
        beware that this is operation has usually no meaning in a series
        >>> p.subs_trunc('x', 3, x=p.lp(1))
        y**2 + 1

        Notice that substitutions are not done one after the other
        >>> p.subs_trunc('x', 3, x=x + y, y=x + 2*y)
        2*x**2 + 6*x*y + 5*y**2
        >>> (x + y)**2 + (x + 2*y)**2
        2*x**2 + 6*x*y + 5*y**2

        which differs from
        >>> p.subs_trunc('x', 3, x=x + y).subs_trunc('x', 3, y=x + 2*y)
        5*x**2 + 12*x*y + 8*y**2
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
            p[zm] = lp.ring_new(n)
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
        """inplace addition of polinomials
        update p1 with values in p2;
        if p1 is a generator, a new polynomial will be returned

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p1 = x + y**2
        >>> p2 = x**2
        >>> p3 = p1
        >>> p1 += p2
        >>> p1
        x**2 + x + y**2
        >>> p1 == p3
        True
        """
        if not p2:
            return p1
        lp1 = p1.lp
        if p1 in lp1.gens:
            return p1 + p2
        if isinstance(p2, LPolyElement):
            if lp1 != p2.lp:
                raise ValueError('p1 and p2 must have the same lp')
            dl = []
            for k, v in p1.iteritems():
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
                p1[mz] = lp1.ring_new(p2)
            else:
                if p2 == -p1[mz]:
                    del p1[mz]
                else:
                    p1[mz] += p2
            return p1

    def __sub__(p1, p2):
        """
        subtract polynomial p2 from p1

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p1 = x + y**2
        >>> p2 = x*y + y**2
        >>> p3 = p1 - p2
        >>> p3
        -x*y + x
        """
        if not p2:
            return p1.copy()
        lp1 = p1.lp
        mz = lp1.zero_mon
        p = LPolyElement(lp1)
        if isinstance(p2, LPolyElement):
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
            elif p2.lp.__class__ == lp1.ring.__class__ and p2.lp == lp1.ring:
                p = p1.copy()
                if mz not in list(p1.keys()):
                    p[mz] = -lp1.ring_new(p2)
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
            p2 = lp1.ring_new(p2)
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

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x + y
        >>> p1 = 4 - p
        >>> p1
        -x - y + 4
        """
        p = LPolyElement(p1.lp)
        for expv in p1:
            p[expv] = -p1[expv]
        p += n
        return p


    def __isub__(p1, p2):
        """
        inplace subtraction of polynomials
        update p1 with values in p2;
        if p1 is a generator, a new polynomial will be returned
        """
        lp1 = p1.lp
        if p1 in lp1.gens:
            return p1 - p2
        if isinstance(p2, LPolyElement):
            if lp1.__class__ != p2.lp.__class__ or lp1 != p2.lp:
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

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p1 = x + y
        >>> p2 = x - y
        >>> p3 = p1*p2
        >>> p3
        x**2 - y**2
        """
        lp1 = p1.lp
        p = LPolyElement(lp1)
        if not p2:
            return p
        if isinstance(p2, LPolyElement):
            if lp1 == p2.lp:
                get = p.get
                p2it = p2.items()
                for exp1, v1 in p1.iteritems():
                    for exp2, v2 in p2it:
                        exp = monomial_mul(exp1, exp2)
                        p[exp] = get(exp, 0) + v1*v2
                p.strip_zero()
                return p
            lp2 = p2.lp
            if lp2.__class__ != lp1.ring.__class__ or lp2 != lp1.ring:
                if lp1.__class__ == lp2.ring.__class__ and lp1 == lp2.ring:
                    p = LPolyElement(p2.lp)
                    for exp2, v2 in p2.iteritems():
                        p[exp2] = p1*v2
                    return p
                else:
                    raise ValueError('p1 and p2 must have the same lp')
        # assume p2 in a number
        for exp1, v1 in p1.iteritems():
            v = v1*p2
            if v:
                p[exp1] = v
        return p

    def __rmul__(p1, p2):
        """
         p2 + p1 with p2 in the coefficient ring of p1

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x + y
        >>> p1 = 4 * p
        >>> p1
        4*x + 4*y
        """
        p = LPolyElement(p1.lp)
        if not isinstance(p2, LPolyElement):
            if not p2:
                return p
        for exp1, v1 in p1.iteritems():
            v = p2*v1
            if v:
                p[exp1] = v
        return p


    def mul_iadd(p, p1, p2):
        """p += p1*p2
        add inplace unless p is a generator, in which case
        return p + p1*p2

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x**2
        >>> p1 = x + y
        >>> p2 = x - y
        >>> p = p.mul_iadd(p1, p2)
        >>> p
        2*x**2 - y**2
        """
        if p in p.lp.gens:
            return p + p1*p2
        if isinstance(p1, LPolyElement) and isinstance(p2, LPolyElement):
            if p1.lp != p2.lp:
                raise ValueError('p1 and p2 must have the same lp')
            get = p.get
            for exp1, v1 in p1.iteritems():
                for exp2, v2 in p2.iteritems():
                    exp = monomial_mul(exp1, exp2)
                    p[exp] = get(exp, 0) + v1*v2
            p.strip_zero()
            return p
        else:
            raise NotImplementedError

    def imul_num(p, c):
        """multiply inplace the polynomial p by an element in the
        coefficient ring, provided p is not one of the generators;
        else multiply not inplace

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x + y**2
        >>> p1 = p.imul_num(3)
        >>> p1
        3*x + 3*y**2
        >>> p1 is p
        True
        >>> p = x
        >>> p1 = p.imul_num(3)
        >>> p1
        3*x
        >>> p1 is p
        False
        """
        if p in p.lp.gens:
            return p*c
        if not c:
            p.clear()
            return
        for exp in p:
            p[exp] *= c
        return p

    def __truediv__(p1, p2):
        """division by a term in the coefficient ring or
        exact division by a polynomial

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p1 = (x**2 + x + y)*(x**2 - y**2)
        >>> p2 = x + y
        >>> p3 = p1/p2
        >>> p4 = (x**2 + x + y)*(x - y)
        >>> p3 == p4
        True
        """
        lp1 = p1.lp
        if isinstance(p2, LPolyElement):
            if lp1 == p2.lp:
                q, r = p1.division([p2])
                if r:
                    raise NotImplementedError('__div__ performs only division without remainder')
                return q[0]
            elif p2.lp.__class__ == lp1.ring.__class__ and p2.lp == lp1.ring:
                zm = p2.lp.zero_mon
                p = LPolyElement(lp1)
                # if p is not a constant, not implemented
                if p2.keys() != [zm]:
                    raise NotImplementedError
                else:
                    p2 = p2[zm]
            else:
                raise NotImplementedError('cannot divide p1 by p2')
        # assume p2 in a number
        p = LPolyElement(lp1)
        if not p2:
            raise ZeroDivisionError
        for exp1, v in p1.iteritems():
            p[exp1] = v/p2
        return p

    __div__ = __truediv__

    def iadd_mon(self, a):
        """add inplace the monomial coeff*x0**i0*x1**i1*...
        a = ((i0, i1, ...), coeff),
        except in the case in which self is a generator, in
        which case add not inplace

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> from sympy.polys.lpoly import monomial_from_sequence
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x**4 + 2*y
        >>> m = monomial_from_sequence((1, 2))
        >>> p1 = p.iadd_mon((m, 5))
        >>> p1
        x**4 + 5*x*y**2 + 2*y
        >>> p1 is p
        True
        >>> p = x
        >>> p1 = p.iadd_mon((m, 5))
        >>> p1
        5*x*y**2 + x
        >>> p1 is p
        False
        """
        if self in self.lp.gens:
            self = self.copy()
        coeff = a[1]
        expv = a[0]
        if expv in self:
            self[expv] += coeff
            if self[expv] == 0:
                del self[expv]
        else:
            self[expv] = coeff
        return self

    def iadd_m_mul_q(p1, p2, mc):
        """ p1 += p2*monom
        m monomial tuple
        c coefficient
        except in the case in which p1 is a generator, in
        which case add not inplace

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> from sympy.polys.lpoly import monomial_from_sequence
        >>> lp, x, y, z = lgens('x, y, z', QQ)
        >>> p1 = x**4 + 2*y
        >>> p2 = y + z
        >>> m = monomial_from_sequence((1, 2, 3))
        >>> p1 = p1.iadd_m_mul_q(p2, (m, 3))
        >>> p1
        x**4 + 3*x*y**3*z**3 + 3*x*y**2*z**4 + 2*y
        """
        if p1 in p1.lp.gens:
            p1 = p1.copy()
        (m, c) = mc
        get = p1.get
        c = p1.lp.ring_new(c)
        for k, v in p2.iteritems():
            ka = monomial_mul(k, m)
            p1[ka] = get(ka, 0) + v*c
        p1.strip_zero()
        return p1

    def leading_expv(self):
        """leading monomial tuple according to the monomial ordering

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y, z = lgens('x, y, z', QQ)
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

        Examples
        ========
       >>> from sympy.polys.domains import QQ
       >>> from sympy.polys.lpoly import lgens
       >>> lp, x, y = lgens('x, y', QQ)
       >>> p = (x + y)**4
       >>> p1 = p.leading_term()
       >>> p1
       x**4
       """
        p = LPolyElement(self.lp)
        expv = self.leading_expv()
        if expv:
            p[expv] = self[expv]
        return p

    def expand(self):
        """expand the coefficients
        """
        if self.lp.SR:
            p =  LPolyElement(self.lp)
            for k, v in self.iteritems():
                p[k] = v.expand()
            return p
        else:
            return self

    def square(p1):
        """square of a polynomial

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x + y**2
        >>> p1 = p.square()
        >>> p1
        x**2 + 2*x*y**2 + y**4
        """
        lp = p1.lp
        if not lp.commuting:
            return p1*p1
        p = LPolyElement(lp)
        get = p.get
        keys = p1.keys()
        for i in range(len(keys)):
            k1 = keys[i]
            pk = p1[k1]
            for j in range(i):
                k2 = keys[j]
                exp = monomial_mul(k1, k2)
                p[exp] = get(exp, 0) + pk*p1[k2]
        p = p.imul_num(2)
        get = p.get
        for k, v in p1.iteritems():
            k2 = monomial_mul(k, k)
            p[k2] = get(k2, 0) + v**2
        p.strip_zero()
        return p

    def pow_miller(p, m):
        """power of a univariate polynomial

        p = sum_i=0**L p_i*x**i
        p**m = sum_k=0**(m*L) a(m, k)*x**k
        Miller pure recurrence formula (see article by
        D. Zeilberger)
        a(m, k) = 1/(k*p_0)*sum_i=1**L p_i*((m + 1)*i - k)*a(m, k - i)

        Reference: D. Zeilberger 'The Miller Recurrence for
        Exponentiating a Polynomial, and its q-Analog'

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p1 = (1 + x + x**2)**5
        >>> p2 = (1 + x + x**2).pow_miller(5)
        >>> p1 == p2
        True
        """

        lp = p.lp
        if lp.ngens > 1:
            raise NotImplementedError("only for univariate polynomials")
        x = lp.gens[0]
        assert lp.ngens == 1
        mindeg = min(p)[0]
        if mindeg != 0:
            p = p/x**mindeg
        degp = max(p)[0]
        pv = [0]*(degp + 1)
        for k, v in p.iteritems():
            pv[k[0]] = v
        a = [lp(0) for i in range(m*degp + 1)]
        a[0] = pv[0]**m
        res = lp(a[0])
        for k in range(1, m*degp + 1):
            s = 0
            for i in range(1, min(degp, k) + 1):
                s += pv[i]*((m + 1)*i - k)*a[k - i]
            a[k] = s/(k*pv[0])
            res[(k, )] = a[k]
        if mindeg != 0:
            res = res*x**(mindeg*m)
        return res

    def __pow__(self, n):
        """raise polynomial to power `n`

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x + y**2
        >>> p1 = p**3
        >>> p1
        x**3 + 3*x**2*y**2 + 3*x*y**4 + y**6
        """
        lp = self.lp
        # test if n is an integer
        n = int_tested(n)
        if n < 0:
            raise ValueError('n >= 0 is required')
        if n == 0:
            if self:
                return lp(1)
            else:
                raise ValueError
        elif len(self) == 1:
            p = LPolyElement(lp)
            k, v = list(self.items())[0]
            # treat case abs(v) = 1 separately to deal with the case
            # in which n is too large to be allowed in v**n
            kn = monomial_pow(k, n)
            if v == 1:
                p[kn] = v
            elif v == -1:
                if n % 2 == 0:
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
        # TODO if lp.SR then use in some cases multinomial coefficients
        if lp.ngens == 1 and n >= 20 and lp.ring in (ZZ, QQ):
            return self.pow_miller(n)
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

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> f = x**3
        >>> f0 = x - y**2
        >>> f1 = x - y
        >>> qv, r = f.division((f0, f1))
        >>> qv[0]
        x**2 + x*y**2 + y**4
        >>> qv[1]
        0
        >>> r
        y**6
        """
        lp = self.lp
        if not self:
            return [], LPolyElement(lp)
        for f in fv:
            if f.lp != lp:
                raise ValueError('self and f must have the same lp')
        s = len(fv)
        qv = [LPolyElement(lp) for i in range(s)]
        p = self.copy()
        r = LPolyElement(lp)
        expvs = [fx.leading_expv() for fx in fv]
        while p:
            i = 0
            divoccurred = 0
            while i < s and divoccurred == 0:
                expv = p.leading_expv()
                expv1 = monomial_div(expv, expvs[i])
                if expv1:
                    c = p[expv]/fv[i][expvs[i]]
                    qv[i] = qv[i].iadd_mon((expv1, c))
                    p = p.iadd_m_mul_q(fv[i], (expv1, -c))
                    divoccurred = 1
                else:
                    i += 1
            if not divoccurred:
                expv =  p.leading_expv()
                r = r.iadd_mon((expv, p[expv]))
                del p[expv]
        if expv == lp.zero_mon:
            r += p
        return qv, r

    def trunc(p1, i, prec):
        """monomials containing x**k, k >= prec neglected
        i is the name of the variable x, or its index

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = (x + y)**4
        >>> p1 = p.trunc('x', 3)
        >>> p1
        6*x**2*y**2 + 4*x*y**3 + y**4
        """
        lp = p1.lp
        if isinstance(i, str):
            i = lp.gens_dict[i]
        p = LPolyElement(lp)
        for exp1 in p1:
            if exp1[i] >= prec:
                continue
            p[exp1] = p1[exp1]
        return p


    def mul_trunc(p1, p2, i, prec):
        """truncation of p1*p2
        p1 and p2 polynomials
        i is the name of the variable x, or its index;
        neglect in p1*p2 the monomials containing x**k, k >= prec

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p1 = (x + y)**2
        >>> p2 = p1.mul_trunc(x**2, x, 3)
        >>> p2
        x**2*y**2
        """
        lp = p1.lp
        if lp.__class__ != p2.lp.__class__ or lp != p2.lp:
            raise ValueError('p1 and p2 must have the same lp')
        if isinstance(i, LPolyElement):
            i = str(i)
        if isinstance(i, str):
            iv = lp.gens_dict[i]
        else:
            iv = i
        p = LPolyElement(p1.lp)
        get = p.get
        items2 = list(p2.items())
        items2.sort(key=lambda e: e[0][iv])
        if lp.ngens == 1:
            for exp1, v1 in p1.iteritems():
                for exp2, v2 in items2:
                    exp = exp1[0] + exp2[0]
                    if exp < prec:
                        exp = (exp, )
                        p[exp] = get(exp, 0) + v1*v2
                    else:
                        break
        else:
            for exp1, v1 in p1.iteritems():
                for exp2, v2 in items2:
                    if exp1[iv] + exp2[iv] < prec:
                        exp = monomial_mul(exp1, exp2)
                        p[exp] = get(exp, 0) + v1*v2
                    else:
                        break
        p.strip_zero()
        return p

    def square_trunc(p1, i, prec):
        """truncation of p1*p1
        i is the name of the variable x, or its index;
        neglect in p1*p1 the monomials containing x**k, k >= prec

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p1 = (x + y)**2
        >>> p2 = p1.square_trunc(x, 3)
        >>> p2
        6*x**2*y**2 + 4*x*y**3 + y**4
        """
        lp = p1.lp
        if not lp.commuting:
            return p1.mul_trunc(p1, i, prec)
        if isinstance(i, LPolyElement):
            i = str(i)
        if isinstance(i, str):
            iv = lp.gens_dict[i]
        else:
            iv = i
        p = LPolyElement(lp)
        get = p.get
        items = list(p1.items())
        items.sort(key=lambda e: e[0][iv])
        for i in range(len(items)):
            exp1, v1 = items[i]
            for j in range(i):
                exp2, v2 = items[j]
                if exp1[iv] + exp2[iv] < prec:
                    exp = monomial_mul(exp1, exp2)
                    p[exp] = get(exp, 0) + v1*v2
                else:
                    break
        p = p.imul_num(2)
        get = p.get
        for expv, v in p1.iteritems():
            if 2*expv[iv] < prec:
                e2 = monomial_mul(expv, expv)
                p[e2] = get(e2, 0) + v**2
        p.strip_zero()
        return p

    def pow_trunc(self, n, i, prec):
        """truncation of self**n, where `n` is rational

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p1 = 1 + x + y
        >>> p2 = p1.pow_trunc(3, y, 2)
        >>> p2
        x**3 + 3*x**2*y + 3*x**2 + 6*x*y + 3*x + 3*y + 1
        """
        lp = self.lp
        i = str(i)
        if isinstance(n, Rational):
            np = n.p
            nq = n.q
            if nq != 1:
                nq = int(sign(np)*nq)
                np =  abs(int(np))
                res = self.nth_root(nq, i, prec)
                if np != 1:
                    res = res.pow_trunc(np, i, prec)
            else:
                res = self.pow_trunc(np, i, prec)
            return res
        assert n == int(n)
        n = int(n)
        if n == 0:
            if self:
                return lp(1)
            else:
                raise ValueError
        if n < 0:
            p1 = self.pow_trunc(-n, i, prec)
            return p1.series_inversion(i, prec)
        if n == 1:
            return self.trunc(i, prec)
        if n == 2:
            return self.square_trunc(i, prec)
        if n == 3:
            p2 = self.square_trunc(i, prec)
            return p2.mul_trunc(self, i, prec)
        p = lp(1)
        while 1:
            if n&1:
                p = self.mul_trunc(p, i, prec)
                n -= 1
                if not n:
                    break
            self = self.square_trunc(i, prec)
            n = n // 2
        return p

    def has_constant_term(p, iv):
        """test if p has a constant term in variable with name iv
        """
        lp = p.lp
        if isinstance(iv, str):
            iv = [lp.gens_dict[iv]]
        else:
            raise NotImplementedError
        zm = lp.zero_mon
        a = [0]*lp.ngens
        for i in iv:
            a[i] = 1
        miv = monomial_from_sequence(a)
        for expv in p:
            if monomial_min(expv, miv) == zm:
                return True
        return False

    def _series_inversion1(p, iv, prec):
        """univariate series inversion 1/p
        iv name of the series variable
        prec precision of the series

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
        for precx in giant_steps(prec):
            tmp = p1.square()
            tmp = tmp.mul_trunc(p, iv, precx)
            p1 = 2*p1 - tmp
        return p1

    def _series_inversion_nc(p, iv, prec):
        """
        a = A**-1, c =  B*a
        (A + B)**-1 = a*(1 - c + c**2 - ...)
        """
        lp = p.lp
        zm = lp.zero_mon
        p0 = p[zm]
        if p0 != lp.ring(1):
          a = p0**-1
        else:
          a = lp.ring(1)
        b = p - p0
        c = -b*a
        cn = lp(1)
        r = lp(1)
        z = lp(0)
        while 1:
          cn = cn.mul_trunc(c, iv, prec)
          if cn == z:
            break
          r += cn
        r = a*r
        return r


    def series_inversion(p, iv, prec):
        """multivariate series inversion 1/p

          iv list of variable names or variable indices
          prec list of truncations for these variables

        In the case of one variable it can also be:

          iv variable name or variable index (0)
          prec truncation integer for the variable

        p is a series with O(x_1**n_1*..x_m**n_m) in
        variables x_k with index or name iv[k - 1]

        p has constant term different from zero

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> from sympy.core import sympify
        >>> from sympy import Symbol
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = (1 + y*x).series_inversion('x', 6)
        >>> p
        -x**5*y**5 + x**4*y**4 - x**3*y**3 + x**2*y**2 - x*y + 1
        >>> a = Symbol('a')
        >>> lp, x = lgens('x', sympify)
        >>> p1 = (1 + a + x).series_inversion('x', 3)
        >>> p1
        ((a + 1)**(-3))*x**2 + (-1/(a + 1)**2)*x + (1/(a + 1))
        """
        lp = p.lp
        zm = lp.zero_mon
        if zm not in p:
            raise NotImplementedError('no constant term in series')
        if (p - p[zm]).has_constant_term(iv):
            raise NotImplementedError('p - p[0] must not have a constant term in the series variables')
        if not lp.commuting:
            return p._series_inversion_nc(iv, prec)
        return p._series_inversion1(iv, prec)

    def derivative(self, n):
        """derivative of p with respect to x_n; the argument n is the
        index of the variable x_n = self.lp.gens[n], or the
        corresponding string.

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x + x**2*y**3
        >>> p1 = p.derivative('x')
        >>> p1
        2*x*y**3 + 1
        """
        lp = self.lp
        if n.__class__ == str:
            n = lp.gens_dict[n]
        else:
            n = int(n)
        p1 = LPolyElement(lp)
        mn = monomial_basis(n, lp.ngens)
        for expv in self:
            if expv[n]:
                e = monomial_div(expv, mn)
                p1[e] = self[expv]*expv[n]
        return p1


    def integrate(self, n):
        """ integrate p with respect to x_n; the argument n is the
        index of the variable x_n = self.lp.gens[n], or the
        corresponding string.

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x + x**2*y**3
        >>> p1 = p.integrate('x')
        >>> p1
        1/3*x**3*y**3 + 1/2*x**2
        """
        lp = self.lp
        if n.__class__ == str:
            n = lp.gens_dict[n]
        p1 = LPolyElement(lp)
        mn = monomial_basis(n, lp.ngens)
        for expv in self:
            e = monomial_mul(expv, mn)
            p1[e] = self[expv]/(expv[n] + 1)
        return p1

############## elementary functions ####################

    def series_from_list(p, c, iv, prec, concur=1):
        """series sum c[n]*p**n
        reduce the number of multiplication summing concurrently
        ax = [1, p, p**2, .., p**(J - 1)]
        s = sum(c[i]*ax[i] for i in range(0, J)) +
            sum(c[i]*ax[i] for i in range(J, 2*J))*p**J +
            sum(c[i]*ax[i] for i in range(2*J, 3*J))*p**(2*J) + ... +
            sum(c[i]*ax[i] for i in range((K - 1)*J, K*J))*p**((K - 1)*J)
        with K >= (n + 1)/J

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> c = [1, 3, 5, 7]
        >>> p1 = (x + y).series_from_list(c, 'x', 3)
        >>> p2 = (1 + 3*(x + y) + 5*(x + y)**2 + 7*(x + y)**3).trunc('x', 3)
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
                q = q.mul_trunc(p, iv, prec)
                s += c[i]*q
            return s
        J = int(math.sqrt(n) + 1)
        if n % J == 0:
            K = n//J
        else:
            K = int(n/J) + 1
        ax = [lp(1)]
        b = 1
        q = lp(1)
        # TODO get a good value
        if len(p) < 20:
            for i in range(1, J):
                q = q.mul_trunc(p, iv, prec)
                ax.append(q)
        else:
            for i in range(1, J):
                if i % 2 == 0:
                    q = ax[i//2].square_trunc(iv, prec)
                else:
                    q = q.mul_trunc(p, iv, prec)
                ax.append(q)
        # optimize using square_trunc
        pj = ax[-1].mul_trunc(p, iv, prec)
        b = lp(1)
        s = lp(0)
        for k in range(K - 1):
            r = J*k
            s1 = c[r]
            for j in range(1, J):
                s1 += c[r + j]*ax[j]
            s1 = s1.mul_trunc(b, iv, prec)
            s += s1
            b = b.mul_trunc(pj, iv, prec)
            if not b:
                break
        k = K - 1
        r = J*k
        if r < n:
            s1 = c[r]*lp(1)
            for j in range(1, J):
                if r + j >= n:
                    break
                s1 += c[r + j]*ax[j]
            s1 = s1.mul_trunc(b, iv, prec)
            s += s1
        return s

    def fun(p, f, *args):
        """
        function of a multivariate series computed by substitution

          p multivariate series
          f method name or function
          args[:-2] arguments of f, apart from the first one
          args[-2] = iv names of the series variables
          args[-1] = prec list of the precisions of the series variables

        The case with f method name is used to compute tan and nth_root
        of a multivariate series:

          p.fun('tan', iv, prec)
          compute _x.tan(iv, prec), then substitute _x with p

          p.fun('nth_root', n, iv, prec)
          compute (_x + p[0]).nth_root(n, '_x', prec)), then substitute _x
          with p - p[0]

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = x + x*y + x**2*y + x**3*y**2
        >>> p.fun('_tan1', 'x', 4)
        1/3*x**3*y**3 + 2*x**3*y**2 + x**3*y + 1/3*x**3 + x**2*y + x*y + x
        """
        lp = p.lp
        lp1 = LPoly(['_x'], lp.ring, lp.order)
        _x = lp1.gens[0]
        h = int(args[-1])
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
        a = sorted(q.iteritems())
        c = [0]*h
        for x in a:
            c[x[0][0]] = x[1]
        p1 = p1.series_from_list(c, args[-2], args[-1])
        return p1

    def _nth_root1(p, n, iv, prec):
        """univariate series nth root of p on commuting ring

          n  integer; compute p**(1/n)
          iv name of the series variable
          prec precision of the series

        The Newton method is used.
        """
        lp = p.lp
        zm = lp.zero_mon
        if zm not in p:
            raise NotImplementedError('no constant term in series')
        n = int_tested(n)
        assert p[zm] == 1
        p1 = lp(1)
        if p == 1:
            return p
        if n == 0:
            return lp(1)
        if n == 1:
            return p
        if n < 0:
            n = -n
            sign = 1
        else:
            sign = 0
        for precx in giant_steps(prec):
            tmp = p1.pow_trunc(n + 1, iv, precx)
            tmp = tmp.mul_trunc(p, iv, precx)
            p1 += p1/n - tmp/n
        if sign:
            return p1
        else:
            return p1._series_inversion1(iv, prec)

    def sqrt(p, iv, prec):
        """square root of multivariate series p

          iv list of variable names or variable indices
          prec list of truncations for these variables

        In the case of one variable it can also be:

          iv variable name or variable index (0)
          prec truncation integer for the variable

        p is a series with O(x_1**n_1*..x_m**n_m) in
        variables x_k with index or name iv[k - 1]

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = (1 + x + x*y).sqrt(x, 3)
        >>> p
        -1/8*x**2*y**2 - 1/4*x**2*y - 1/8*x**2 + 1/2*x*y + 1/2*x + 1
        """
        iv = str(iv)
        p1 = p.nth_root(-2, iv, prec)
        return p.mul_trunc(p1, iv, prec)


    def nth_root(p, n, iv, prec):
        """multivariate series nth root of p

          n  integer; compute p**(1/n)
          iv list of variable names or variable indices
          prec list of truncations for these variables

        In the case of one variable it can also be:

          iv variable name or variable index (0)
          prec truncation integer for the variable

        p is a series with O(x_1**n_1*..x_m**n_m) in
        variables x_k with index or name iv[k - 1]

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x, y = lgens('x, y', QQ)
        >>> p = (1 + x + x*y).nth_root(-3, x, 3)
        >>> p
        2/9*x**2*y**2 + 4/9*x**2*y + 2/9*x**2 - 1/3*x*y - 1/3*x + 1
        """
        if n == 0 and p == 0:
                raise ValueError('0**0 expression')
        iv = str(iv)
        if n == 1:
            return p.trunc(iv, prec)
        lp = p.lp
        if (p - 1).has_constant_term(iv):
            if lp.zero_mon in p:
                c = p[lp.zero_mon]
                if isinstance(c, PythonRationalType):
                    c1 = Rational(c.p, c.q)
                    cn = Pow(c1, S.One/n)
                else:
                    cn = Pow(c, S.One/n)
                if cn.is_Rational:
                    if not lp.SR:
                        cn = lp.ring(cn.p, cn.q)
                    return cn*(p/c).nth_root(n, iv, prec)

            if not lp.SR:
                raise TaylorEvalError('p - 1 must not have a constant term in the series variables')
            else:
                if lp.zero_mon in p:
                    c = p[lp.zero_mon]
                    if c.is_positive:
                        return (p/c).nth_root(n, iv, prec)*c**Rational(1, n)
                    else:
                        raise NotImplementedError
        if lp.commuting and lp.ngens == 1:
            return p._nth_root1(n, iv, prec)
        else:
            return p.fun('_nth_root1',n, iv, prec)

    def re_acoth(p, iv, prec):
        """Re(acoth(p)) = 1/2*(log(1 + p) - log(1 - p))

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = x.re_acoth(x, 8)
        >>> p
        1/7*x**7 + 1/5*x**5 + 1/3*x**3 + x
        """
        iv = str(iv)
        if p.has_constant_term(iv):
            raise NotImplementedError('p must not have a constant term in the series variables')
        p = p
        lp = p.lp
        s = p.copy()
        q = p
        n = 3
        z = lp(0)
        p2 = p.square_trunc(iv, prec)
        while 1:
            q = q.mul_trunc(p2, iv, prec)
            if q == z:
                break
            s += q/n
            n += 2
        return s

    def acot1(p, iv, prec):
        """acot1(p) = acot(p) - constant part = -p +p**3/3 -p**5/5 + ...
        """
        iv = str(iv)
        if p.has_constant_term(iv):
            raise NotImplementedError('p must not have a constant term in the series variables')
        p = p
        lp = p.lp
        s = -p
        q = -p
        n = 3
        z = lp(0)
        p2 = -p.square_trunc(iv, prec)
        while 1:
            q = q.mul_trunc(p2, iv, prec)
            if q == z:
                break
            s += q/n
            n += 2
        return s

    def _log_series(p, iv, prec):
        lp = p.lp
        one = lp.ring(1)
        c = [one*0]
        for k in range(1, prec):
            if k % 2 == 0:
                cf = -one/k
            else:
                cf = one/k
            c.append(cf)
        return (p - 1).series_from_list(c, iv, prec)

    def log(p, iv, prec):
        """
        logarithm of p with truncation

        For univariate series or with the total degree
        truncation integral dx p**-1*d p/dx is used.

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = (1 + x).log(x, 8)
        >>> p
        1/7*x**7 - 1/6*x**6 + 1/5*x**5 - 1/4*x**4 + 1/3*x**3 - 1/2*x**2 + x
        """
        lp = p.lp
        iv = str(iv)
        if (p - 1).has_constant_term(iv):
            if not lp.SR:
                raise TaylorEvalError('p - 1 must not have a constant term in the series variables')
            else:
                if lp.zero_mon in p:
                    c = p[lp.zero_mon]
                    if c.is_positive:
                        return sympy.functions.log(c) + (p/c).log(iv, prec)
                    else:
                        raise NotImplementedError
            raise NotImplementedError('p - 1 must not have a constant term in the series variables')
        lp = p.lp
        if lp.commuting:
            dlog = p.derivative(iv)
            dlog = dlog.mul_trunc(p._series_inversion1(iv, prec), iv, prec - 1)
            return dlog.integrate(iv)
        return p._log_series(iv, prec)

    def _atan_series(p, iv, prec):
        lp = p.lp
        mo = lp.ring(-1)
        c = [-mo]
        p2 = p.square_trunc(iv, prec)
        for k in range(1, prec):
            c.append(mo**k/(2*k + 1))
        s = p2.series_from_list(c, iv, prec)
        s = s.mul_trunc(p, iv, prec)
        return s


    def atan(p, iv, prec):
        """
        arctangent of a series

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = x.atan(x, 8)
        >>> p
        -1/7*x**7 + 1/5*x**5 - 1/3*x**3 + x
        """
        iv = str(iv)
        if p.has_constant_term(iv):
            raise NotImplementedError('polynomial must not have constant term in the series variables')
        lp = p.lp
        if iv in lp.pol_gens and lp.commuting:
            dp = p.derivative(iv)
            p1 = p.square_trunc(iv, prec) + 1
            p1 = p1.series_inversion(iv, prec - 1)
            p1 = dp.mul_trunc(p1, iv, prec - 1)
            return p1.integrate(iv)
        else:
            return p._atan_series(iv, prec)

    def lambert(p, iv, prec):
        """
        principal branch of the Lambert W function

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = x.lambert(x, 8)
        >>> p
        16807/720*x**7 - 54/5*x**6 + 125/24*x**5 - 8/3*x**4 + 3/2*x**3 - x**2 + x
        """
        lp = p.lp
        iv = str(iv)
        p1 = lp(0)
        if p.has_constant_term(iv):
            raise NotImplementedError('polynomial must not have constant term in the series variables')
        if iv in lp.pol_gens and lp.commuting:
            for precx in giant_steps(prec):
                e = p1.exp(iv, precx)
                p2 = e.mul_trunc(p1, iv, precx) - p
                p3 = e.mul_trunc(p1 + 1, iv, precx)
                p3 = p3.series_inversion(iv, precx)
                tmp = p2.mul_trunc(p3, iv, precx)
                p1 -= tmp
            return p1
        else:
            raise NotImplementedError

    def asin(p, iv, prec):
        """
        arcsin of a series

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = x.asin(x, 8)
        >>> p
        5/112*x**7 + 3/40*x**5 + 1/6*x**3 + x
        """
        iv = str(iv)
        if p.has_constant_term(iv):
            raise NotImplementedError('polynomial must not have constant term in the series variables')
        lp = p.lp
        if iv in lp.pol_gens and lp.commuting:
            # get a good value
            if len(p) > 20:
                dp = p.derivative(iv)
                p1 = 1 - p.square_trunc(iv, prec - 1)
                p1 = p1.nth_root(-2, iv, prec - 1)
                p1 = dp.mul_trunc(p1, iv, prec - 1)
                return p1.integrate(iv)
            one = lp.ring(1)
            c = [0, one, 0]
            for k in range(3, prec, 2):
                if k % 2 == 1:
                    c.append((k - 2)**2*c[-2]/(k*(k - 1)))
                    c.append(0)
            return p.series_from_list(c, iv, prec)

        else:
            raise NotImplementedError

    def asinh(p, iv, prec):
        """
        arcsinh of a series

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = x.asinh(x, 8)
        >>> p
        -5/112*x**7 + 3/40*x**5 - 1/6*x**3 + x
        """
        iv = str(iv)
        if p.has_constant_term(iv):
            raise NotImplementedError('polynomial must not have constant term in the series variables')
        lp = p.lp
        if iv in lp.pol_gens and lp.commuting:
            dp = p.derivative(iv)
            p1 = 1 + p.square_trunc(iv, prec - 1)
            p1 = p1.nth_root(-2, iv, prec - 1)
            p1 = dp.mul_trunc(p1, iv, prec - 1)
            return p1.integrate(iv)
        else:
            raise NotImplementedError

    def _atanh_series(p, iv, prec):
        lp = p.lp
        one = lp.ring(1)
        c = [one]
        p2 = p.square_trunc(iv, prec)
        for k in range(1, prec):
            c.append(one/(2*k + 1))
        s = p2.series_from_list(c, iv, prec)
        s = s.mul_trunc(p, iv, prec)
        return s

    def atanh(p, iv, prec):
        """ hyperbolic arctangent

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = x.atanh(x, 8)
        >>> p
        1/7*x**7 + 1/5*x**5 + 1/3*x**3 + x
        """
        iv = str(iv)
        if p.has_constant_term(iv):
            raise NotImplementedError('polynomial must not have constant term in the series variables')
        lp = p.lp
        if iv in lp.pol_gens and lp.commuting:
            dp = p.derivative(iv)
            p1 = -p.square_trunc(iv, prec) + 1
            p1 = p1.series_inversion(iv, prec - 1)
            p1 = dp.mul_trunc(p1, iv, prec - 1)
            return p1.integrate(iv)
        else:
            return p._atanh_series(iv, prec)

    def _tanh1(p, iv, prec):
        lp = p.lp
        p1 = lp(0)
        for precx in giant_steps(prec):
            tmp = p - p1.atanh(iv, precx)
            tmp = tmp.mul_trunc(1 - p1.square(), iv, precx)
            p1 += tmp
        return p1

    def tanh(p, iv, prec):
        """ hyperbolic tangent of a series

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = x.tanh(x, 8)
        >>> p
        -17/315*x**7 + 2/15*x**5 - 1/3*x**3 + x
        """
        lp = p.lp
        iv = str(iv)
        if p.has_constant_term(iv):
            raise NotImplementedError('p must not have constant part in series variables')
        if lp.commuting and lp.ngens == 1:
            return p._tanh1(iv, prec)
        return p.fun('_tanh1', iv, prec)

    def _tan1(p, iv, prec):
        lp = p.lp
        p1 = lp(0)
        for precx in giant_steps(prec):
            tmp = p - p1.atan(iv, precx)
            tmp = tmp.mul_trunc(1 + p1.square(), iv, precx)
            p1 += tmp
        return p1

    def tan(p, iv, prec):
        """tangent of a series

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = x.tan(x, 8)
        >>> p
        17/315*x**7 + 2/15*x**5 + 1/3*x**3 + x
        """
        iv = str(iv)
        lp = p.lp
        if p.has_constant_term(iv):
            raise NotImplementedError('p must not have constant part in series variables')
        if lp.commuting and lp.ngens == 1:
            return p._tan1(iv, prec)
        return p.fun('tan', iv, prec)


    def _exp1(p, iv, prec):
        """
        exponential of a univariate series, or of a multivariate
        series with total degree truncation
        """
        lp = p.lp
        zm = lp.zero_mon
        p1 = lp(1)
        for precx in giant_steps(prec):
            tmp = (p - p1.log(iv, precx)).mul_trunc(p1, iv, precx)
            p1 += tmp
        return p1

    def exp(p, iv, prec):
        """
        exponential of a series

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = x.exp(x, 6)
        >>> p
        1/120*x**5 + 1/24*x**4 + 1/6*x**3 + 1/2*x**2 + x + 1
        """
        lp = p.lp
        iv = str(iv)
        if p.has_constant_term(iv):
            zm = lp.zero_mon
            if not lp.SR:
                raise TaylorEvalError('p must not have constant part in series variables')
            return sympy.functions.exp(p[zm])*(p - p[zm]).exp(iv, prec)
        if len(p) > 20 and lp.commuting:
            return p._exp1(iv, prec)
        one = lp.ring(1)
        n = 1
        k = 1
        c = []
        for k in range(prec):
            c.append(one/n)
            k += 1
            n *= k

        r = p.series_from_list(c, iv, prec)
        return r

    def sin(p, iv, prec):
        """ sin of a series

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = x.sin(x, 6)
        >>> p
        1/120*x**5 - 1/6*x**3 + x
        """
        lp = p.lp
        iv = str(iv)
        if p.has_constant_term(iv):
            if not lp.SR:
                raise TaylorEvalError
            zm = lp.zero_mon
            c = p[zm]
            if c.is_number and not c.is_real:
                raise TaylorEvalError
            p1 = p - c
            return sympy.functions.cos(c)*p1.sin(iv, prec) + sympy.functions.sin(c)*p1.cos(iv, prec)
        # get a good value
        if len(p) > 20 and lp.ngens == 1:
            t = (p/2).tan(iv, prec)
            t2 = t.square_trunc(iv, prec)
            p1 = (1 + t2).series_inversion(iv, prec)
            return p1.mul_trunc(2*t, iv, prec)
        one = lp.ring(1)
        n = 1
        c = [0]
        for k in range(2, prec + 2, 2):
            c.append(one/n)
            c.append(0)
            n *= -k*(k + 1)
        return p.series_from_list(c, iv, prec)

    def cos(p, iv, prec):
        """ cos of a series

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = x.cos(x, 6)
        >>> p
        1/24*x**4 - 1/2*x**2 + 1
        """
        lp = p.lp
        iv = str(iv)
        if p.has_constant_term(iv):
            zm = lp.zero_mon
            if not lp.SR:
                raise TaylorEvalError
            c = p[zm]
            if not c.is_real:
                raise NotImplementedError
            p1 = p - c
            return sympy.functions.cos(c)*p1.cos(iv, prec) - \
                    sympy.functions.sin(c)*p1.sin(iv, prec)
        # get a good value
        if len(p) > 20 and lp.ngens == 1:
            t = (p/2).tan(iv, prec)
            t2 = t.square_trunc(iv, prec)
            p1 = (1 + t2).series_inversion(iv, prec)
            return p1.mul_trunc(1 - t2, iv, prec)
        one = lp.ring(1)
        n = 1
        c = []
        for k in range(2, prec + 2, 2):
            c.append(one/n)
            c.append(0)
            n *= -k*(k - 1)
        return p.series_from_list(c, iv, prec)


    def cos_sin(p, iv, prec):
        """
        tuple (p.cos(iv, iv), p.sin(iv, iv))
        """
        iv = str(iv)
        t = (p/2).tan(iv, prec)
        t2 = t.square_trunc(iv, prec)
        p1 = (1 + t2).series_inversion(iv, prec)
        return (p1.mul_trunc(1 - t2, iv, prec), p1.mul_trunc(2*t, iv, prec))

    def sinh(self, iv, prec):
        """ hyperbolic sin of a series

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = x.sinh(x, 8)
        >>> p
        1/5040*x**7 + 1/120*x**5 + 1/6*x**3 + x
        """
        iv = str(iv)
        t = self.exp(iv, prec)
        t1 = t.series_inversion(iv, prec)
        return (t - t1)/2

    def cosh(p, iv, prec):
        """ cos of a series

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> lp, x = lgens('x', QQ)
        >>> p = x.cosh(x, 8)
        >>> p
        1/720*x**6 + 1/24*x**4 + 1/2*x**2 + 1
        """
        iv = str(iv)
        t = p.exp(iv, prec)
        t1 = t.series_inversion(iv, prec)
        return (t + t1)/2

    def cosh_sinh(p, iv, prec):
        """ tuple (p.cosh(iv, iv), p.sinh(iv, iv))
        """
        iv = str(iv)
        t = p.exp(iv, prec)
        t1 = t.series_inversion(iv, prec)
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
        """Convert a LPolyElement into a Sympy expression.

        Examples
        ========
        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.lpoly import lgens
        >>> from sympy import Symbol
        >>> lp, X = lgens('X', QQ)
        >>> x = Symbol('x')
        >>> p = (1 + X)**3
        >>> p
        X**3 + 3*X**2 + 3*X + 1
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
        for monom, coeff in p.iteritems():
            term = [convert(coeff)]
            for g, m in zip(gens, monom):
                term.append(Pow(g, m))

            result.append(Mul(*term))
        return Add(*result)

class LPolySubs(object):
    """class for substitutions of variables with polynomials,
    possibly truncated.

    Examples
    ========
    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens
    >>> from sympy.polys.lpoly import LPolySubs
    >>> lp, x, y = lgens('x, y', QQ)
    >>> sb = LPolySubs(lp, lp, {'y':y + 1})
    >>> p1 = sb.subs((x + y)**2)
    >>> p1
    x**2 + 2*x*y + 2*x + y**2 + 2*y + 1
    >>> p1 = sb.subs_trunc((x + y)**4, 'x', 2)
    >>> p1
    4*x*y**3 + 12*x*y**2 + 12*x*y + 4*x + y**4 + 4*y**3 + 6*y**2 + 4*y + 1
    """
    def __init__(self, lp1, lp2, rules):
        """
        initialize a substitution object

        lp1, lp2 LPoly objects

        rules    rules of substitution x=px, y=py, .. where x, y are
                 in p.lp.gens and px, py are the polynomials with which
                 they are substituted:
                 p = ... + c_{i,j}*x**i*y**j + ... ->
                     ... + c_{i,j}*px**i*py**j + ...

        TODO currently it must be lp1.ring == lp2.ring
             Implement substitutions between different polynomial rings,
             e.g. mpq and mpf
        """
        self.lp1 = lp1
        gens_dict = lp1.gens_dict
        self.lp2 = lp2
        if lp1.ring.__class__ != lp2.ring.__class__ or lp1.ring != lp2.ring:
            raise NotImplementedError
        d = {} # replace monomials with (i, pw)
        gens = lp1.gens
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
        p1 = LPolyElement(lp2)
        pk = sorted(p.keys())
        for expv in pk:
            p2 = lp2(1)
            for i in range(ngens):
                pw = expv[i]
                if pw == 0:
                    continue
                if (i, pw) not in d:
                    w, r = divmod(pw, 2)
                    if r == 0 and (i, w) in d:
                        d[(i, pw)] = d[(i, w)]**2
                    elif (i, pw - 1) in d:
                        d[(i, pw)] = d[(i, pw - 1)]*d[(i, 1)]
                    else:
                        d[(i, pw)] = d[(i, 1)]**pw
                p2 *= d[(i, pw)]
            p1 += p2*p[expv]
        return p1

    def subs_trunc(self, p, ii, h):
        """substitution with truncation in variable ii with precision h

          p     input polynomial
          iv    (name of the) variable in which the series truncation is done
          nv    order of the truncation
        """
        ii = str(ii)
        lp1 = self.lp1
        lp2 = self.lp2
        ngens = lp1.ngens
        assert p.lp == self.lp1
        d = self.d.copy()
        p1 = lp2(0)
        pk = sorted(p.keys())
        for expv in pk:
            p2 = lp2(1)
            for i in range(ngens):
                pw = expv[i]
                if pw == 0:
                    continue
                if (i, pw) not in d:
                    w, r = divmod(pw, 2)
                    if r == 0 and (i, w) in d:
                        d[(i, pw)] = d[(i, w)].square_trunc(ii, h)
                    elif (i, pw - 1) in d:
                        d[(i, pw)] = d[(i, pw - 1)].mul_trunc(d[(i, 1)], ii, h)
                    else:
                        d[(i, pw)] = d[(i, 1)].pow_trunc(pw, ii, h)
                p2 = p2.mul_trunc(d[(i, pw)], ii, h)
            p1 += p2*p[expv]
        return p1
