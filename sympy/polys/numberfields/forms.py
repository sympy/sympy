"""Formal representations of algebraic number theoretic objects. """

from sympy import Rational, igcd, ilcm, sympify
from sympy.polys import Poly
from sympy.polys.domains import GF, QQ, ZZ
from sympy.polys.numberfields.modules import (
    HnfEndomorphismRing, ModuleElement, ModuleWithDenominator, Submodule)
from sympy.polys.numberfields.utilities import is_int, is_rat, is_zero
from sympy.utilities import public


@public
class StandardRep:
    """
    The "standard representation" for an algebraic number.

    Consists of:
      - a monic irreducible polynomial T(x) over ZZ
      - a tuple (c0, c1, ..., cn-1) of integers, where n = deg(T)
      - a positive integer d
    and represents the number
        sum(c[i]*theta**i/d for i in range(n))
    where theta is a root of T.
    """

    def __init__(self, T, coeffs, denom):
        """
        Parameters
        ----------
        T: monic, irreducible polynomial over ZZ defining the number field
          K(theta) to which this algebraic number belongs, T(theta) = 0.
        coeffs: iterable of length n == deg(T), giving the integer coefficients
          of the powers 1, theta, ..., theta**(n-1) in this number's numerator.
        denom: positive integer giving the denominator of this number.

        Note: To be in standard form, should be gcd(*coeffs, denom) = 1.
        This class offers methods to check this condition, and to create this
        condition if necessary, but it is up to you to call these methods.

        The arithmetic methods (__add__, __mul__, etc.) are designed to maintain
        this condition, provided it is present to begin with.
        """
        self.T = T
        self.coeffs = list(coeffs)
        self.denom = denom
        self.n = len(self.coeffs)
        self._num_gcd = None
        self._inverse = None

    @property
    def num_gcd(self):
        """Check the GCD of the numerator coefficients."""
        if self._num_gcd is None:
            if self.n == 1:
                self._num_gcd = abs(self.coeffs[0])
            else:
                self._num_gcd = igcd(*self.coeffs)
        return self._num_gcd

    def is_reduced(self):
        g = self.num_gcd
        return igcd(g, self.denom) == 1

    def reduced(self):
        """Produce a reduced version of this rep."""
        g = igcd(self.denom, *self.coeffs)
        return StandardRep(self.T, [c // g for c in self.coeffs], self.denom // g)

    def clone(self):
        return StandardRep(self.T, self.coeffs[:], self.denom)

    standard_rep = clone

    def in_rationals(self):
        return all(c == 0 for c in self.coeffs[1:])

    @property
    def is_zero(self):
        return self == 0

    @classmethod
    def zero(cls, T):
        return cls.from_int(T, 0)

    @classmethod
    def one(cls, T):
        return cls.from_int(T, 1)

    @classmethod
    def from_int(cls, T, a):
        return cls(T, [a] + [0] * (T.degree() - 1), 1)

    @classmethod
    def from_rational(cls, T, a):
        a = sympify(a)
        n, d = a.as_numer_denom()
        return cls(T, [n] + [0] * (T.degree() - 1), d)

    @classmethod
    def from_poly(cls, T, f):
        n, k = T.degree(), f.degree()
        if k >= n:
            f = f % T
        if f == 0:
            return cls.zero(T)
        d, g = f.clear_denoms()
        c = list(reversed(g.all_coeffs()))
        ell = len(c)
        z = [ZZ(0)] * (n - ell)
        return cls(T, c + z, d)

    @classmethod
    def from_matrix_col(cls, T, M, j, d):
        return cls(T, M.col(j).flat(), d)

    def __eq__(self, other):
        if isinstance(other, StandardRep):
            s_r = self if self.is_reduced() else self.reduced()
            o_r = other if other.is_reduced() else other.reduced()
            return s_r.denom == o_r.denom and s_r.coeffs == o_r.coeffs
        if is_rat(other):
            other = sympify(other)
            return self.in_rationals() and other.numerator * self.denom == other.denominator * self.coeffs[0]
        return NotImplemented

    def __add__(self, other):
        if isinstance(other, StandardRep):
            d, e = self.denom, other.denom
            m = ilcm(d, e)
            u, v = m // d, m // e
            r = StandardRep(
                self.T, [u * a + v * b for a, b in zip(self.coeffs, other.coeffs)], m)
            return r.reduced()
        if is_rat(other):
            other = sympify(other)
            return self + StandardRep.from_rational(self.T, other)
        return NotImplemented

    __radd__ = __add__

    def __neg__(self):
        return self * (-1)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, r):
        if isinstance(r, StandardRep):
            f, g = self.poly(), r.poly()
            e, h = (f * g % self.T).clear_denoms()
            return StandardRep.from_poly(self.T, h) // e
        if is_rat(r):
            r = sympify(r)
            a, b = r.as_numer_denom()
            s = StandardRep(self.T, [a*c for c in self.coeffs], b*self.denom)
            return s.reduced()
        return NotImplemented

    __rmul__ = __mul__

    def __pow__(self, a):
        if is_int(a):
            if a < 0:
                return self.inverse ** (-a)
            elif a == 0:
                return StandardRep.one(self.T)
            elif a == 1:
                return self.clone()
            else:
                f = self.numerator()
                # TODO: check whether this is an efficient pow-mod operation or not:
                g = f ** a % self.T
                return StandardRep.from_poly(self.T, g) // (self.denom ** a)
        return NotImplemented

    def __floordiv__(self, a):
        if is_rat(a):
            a = sympify(a)
            if a == 0:
                raise ZeroDivisionError
            return self * (1/a)
        return NotImplemented

    @property
    def inverse(self):
        if self._inverse is None:
            f = self.poly()
            e, h = f.invert(self.T).clear_denoms()
            self._inverse = StandardRep.from_poly(self.T, h) // e
        return self._inverse

    def __rfloordiv__(self, a):
        return a * self.inverse

    def __mod__(self, a):
        """
        Reducing mod an integer a reduces all numerator coeffs mod d*a, where
        d is our denominator.

        For example, if we represent

            15*x + 1     x + 1
            --------  =  ----- + 7*x
               2           2

        then reducing mod 7 should mean throwing away that part that is a poly
        divisible by 7. But this is achieved by reducing our coeffs mod 2*7.
        """
        if is_int(a):
            a = sympify(a)
            m = a * self.denom
            return StandardRep(self.T, [b % m for b in self.coeffs], self.denom)
        return NotImplemented

    def numerator(self, x=None):
        """Obtain the numerator as a polynomial over ZZ."""
        x = x or self.T.gen
        return Poly(reversed(self.coeffs), x, domain=ZZ)

    def poly(self, x=None):
        """Obtain the number as a polynomial over QQ."""
        return self.numerator(x=x) // self.denom

    def norm(self, T=None):
        """Compute the norm of this number."""
        T = T or self.T
        x = T.gen
        A = self.numerator(x=x)
        return T.resultant(A) // self.denom**self.n


@public
class HNF(ModuleWithDenominator):
    """
    Hermite Normal Form representation for orders and ideals.

    Consists of:
      - a monic irreducible polynomial T(x) over ZZ
      - a square, upper triangular n x n matrix W with integer entries,
        where n = deg(T)
      - a positive integer d
    Each column equates to the standard representation of an algebraic number,
    with the single denominator d being common to all columns. The representations
    are thus over the power basis in theta, a root of T.
    """

    def __init__(self, T, W, d):
        if not W.is_square or not W.is_upper:
            raise ValueError('HNF matrix must be square and upper triangular')
        super().__init__(W, d)
        self.T = T
        self._mult_tab = None

    @classmethod
    def for_power_basis(cls, T):
        from sympy.matrices import Matrix
        return cls(T, Matrix.eye(T.degree()), 1)

    def __str__(self):
        return str((self.W, self.d))

    def __eq__(self, other):
        if isinstance(other, HNF):
            return self.QQ_matrix() == other.QQ_matrix()
        return NotImplemented

    @property
    def mult_tab(self):
        if self._mult_tab is None:
            omegas = self.standard_reps()
            n = self.n
            M = {}
            for u in range(n):
                M[u] = {}
                om_u = omegas[u]
                for v in range(u, n):
                    om_v = omegas[v]
                    c = self.represent(om_u * om_v, domain=ZZ)
                    M[u][v] = c
            self._mult_tab = M
        return self._mult_tab

    @classmethod
    def zero(cls, T):
        return cls.from_rational(T, 0)

    @classmethod
    def one(cls, T):
        return cls.from_rational(T, 1)

    @classmethod
    def from_rational(cls, T, a):
        n = T.degree()
        from sympy.matrices import Matrix
        return cls(T, a * Matrix.eye(n), 1)

    def __add__(self, other):
        if isinstance(other, ModuleWithDenominator):
            d, e = self.d, other.d
            m = ilcm(d, e)
            a, b = m // d, m // e
            W = (a * self.W).row_join(b * other.W)
            from sympy.matrices.normalforms import hermite_normal_form
            W = hermite_normal_form(W)
            return HNF(self.T, W, m)
        return NotImplemented

    __radd__ = __add__

    def is_compatible_module_element(self, c):
        if isinstance(c, ModuleElement):
            if isinstance(c.module, Submodule):
                return c.module.container == self
        return False

    def __mul__(self, c):
        if is_zero(c) or (isinstance(c, StandardRep) and c == 0):
            return HNF.zero(self.T)
        elif is_rat(c):
            return self * StandardRep.from_rational(self.T, c)
        elif isinstance(c, StandardRep):
            omegas = self.standard_reps()
            cols = [c * om for om in omegas]
            return HNF.from_std_reps(cols)
        elif isinstance(c, HNF):
            alphas, betas = self.standard_reps(), c.standard_reps()
            cols = [a * b for a in alphas for b in betas]
            return HNF.from_std_reps(cols)
        elif self.is_compatible_module_element(c):
            return NotImplemented  # TODO
        return NotImplemented

    __rmul__ = __mul__

    def det(self):
        D = 1
        for i in range(self.n):
            D *= self.W[i, i]
        return D // self.d**self.n

    def QQ_matrix(self):
        return self.W / self.d

    @classmethod
    def from_std_reps(cls, R):
        """
        Build an HNF from a list of StandardReps.
        All reps must be a common length n, and there must be n or more of them.
        The reps need not already form a basis; the Hermite normal form reduction
        algorithm will be applied in order to try to form a basis.
        """
        from sympy.matrices import Matrix
        n = len(R)
        if not n > 0:
            raise ValueError('Need at least one number!')
        m = R[0].n
        if not all(r.n == m for r in R):
            raise ValueError('Standard reps must have same length for HNF.')
        from sympy.matrices.normalforms import hermite_normal_form
        d = R[0].denom if n == 1 else ilcm(*[r.denom for r in R])
        A = [d // r.denom for r in R]
        M = Matrix([[a*c for c in r.coeffs] for a, r in zip(A,R)]).T
        W = hermite_normal_form(M)
        return cls(R[0].T, W, d)

    def represent(self, r, domain=None, modulus=None):
        """
        Represent a given StandardRep as a linear combination over this basis.

        Parameters
        ----------
        r: StandardRep instance. Length must equal number of cols in this HNF.
        domain: (optional) QQ or ZZ to indicate that the coefficients of the
          linear combination must belong to this domain.
        modulus: (optional) a positive rational prime p to indicate that the
          coefficients of the linear combination must belong to GF(p).

        If neither domain nor modulus is specified then we assume that the
        coeffs may belong to QQ.

        If ZZ or GF(p) is specified and any coefficient fails to belong to this
        domain, we raise `ValueError`.

        Returns
        -------
        List of coefficients, belonging to QQ, ZZ, or GF(p), as specified above.
        The jth coefficient corresponds to the jth column of this HNF.

        Raises
        ------
        ValueError if dimensions are mismatched or if the given number is not
        representable over this basis, and with coeffs in the desired domain.
        """
        # Coerce ModuleElement, if need be:
        r = self.standard_rep(r)
        if r.n != self.n:
            raise ValueError('Mismatched HNF and StandardRep dimensions!')
        R = QQ
        if domain is not None:
            if domain not in [QQ, ZZ]:
                raise ValueError('Domain must be QQ or ZZ.')
            R = domain
        elif modulus is not None:
            if modulus < 2:  # skip primality test for sake of speed
                raise ValueError('If defined, modulus must be positive prime.')
            R = GF(modulus)
        n = self.n
        target = r
        coeffs = [R(0)] * n
        for i in range(1, n + 1):
            a = target.coeffs[n - i]
            if a != 0:
                col = self.standard_rep(n - i, reduced=False)
                quo = Rational(a, target.denom) / Rational(col.coeffs[n - i], col.denom)
                if R.is_ZZ:
                    if not is_int(quo):
                        raise ValueError('Outside ZZ-span of basis!')
                    coeffs[n - i] = R(quo)
                elif R.is_FF:
                    numer, denom = quo.as_numer_denom()
                    if R(denom) == 0:
                        raise ValueError(f'Outside {R}-span of basis!')
                    coeffs[n - i] = R(numer)/R(denom)
                else:
                    coeffs[n - i] = R(quo)
                target = target - quo * col
        assert target == 0
        return coeffs

    def standard_rep(self, j, reduced=False):
        if isinstance(j, StandardRep):
            return j
        if isinstance(j, ModuleElement):
            col = self.W * j.column()
        else:
            col = self.W.col(j)
        r = StandardRep(self.T, col.flat(), self.d)
        if reduced:
            r = r.reduced()
        return r

    def standard_reps(self, reduced=False):
        return [self.standard_rep(j, reduced=reduced) for j in range(self.n)]

    def endomorphism_ring(self):
        return HnfEndomorphismRing(self)
