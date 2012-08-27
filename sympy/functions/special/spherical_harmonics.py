from sympy import pi, I
from sympy.core.basic import C
from sympy.core.singleton import S
from sympy.core import Dummy, sympify
from sympy.core.function import Function, ArgumentIndexError
from sympy.functions import legendre, assoc_legendre
from sympy.functions.elementary.miscellaneous import sqrt

Pl = legendre
Plm = assoc_legendre

_x = Dummy("x")

class Ynm(Function):
    r"""
    Spherical harmonics defined as

    .. math::
        Y_n^m(\theta, \varphi) := \sqrt{\frac{(2n+1)(n-m)!}{4\pi(n+m)}}
                                  \exp(-i m \varphi)
                                  \mathrm{P}_n^m\left(\cos(\theta)\right)

    Ynm() gives the spherical harmonic function of order `n` and `m`
    in `\theta` and `\varphi`, `Y_n^m(\theta, \varphi)`. The four
    parameters are as follows: `n \geq 0` an integer and `m` an integer
    such that `-n \leq m \leq n` holds. The two angles are real-valued
    with `\theta \in [0, \pi]` and `\varphi \in [0, 2\pi]`.

    Examples
    ========

    >>> from sympy import Ynm, Symbol
    >>> from sympy.abc import n,m
    >>> theta = Symbol("theta")
    >>> phi = Symbol("phi")

    >>> Ynm(n, m, theta, phi)
    Ynm(n, m, theta, phi)

    Several symmetries are known, for the order

    >>> from sympy import Ynm, Symbol
    >>> from sympy.abc import n,m
    >>> theta = Symbol("theta")
    >>> phi = Symbol("phi")

    >>> Ynm(n, -m, theta, phi)
    (-1)**m*exp(-2*I*m*phi)*Ynm(n, m, theta, phi)

    as well as for the angles

    >>> from sympy import Ynm, Symbol, simplify
    >>> from sympy.abc import n,m
    >>> theta = Symbol("theta")
    >>> phi = Symbol("phi")

    >>> Ynm(n, m, -theta, phi)
    Ynm(n, m, theta, phi)

    >>> Ynm(n, m, theta, -phi)
    exp(-2*I*m*phi)*Ynm(n, m, theta, phi)

    For specific integers n and m we can evalute the harmonics
    to more useful expressions

    >>> simplify(Ynm(0, 0, theta, phi).expand(func=True))
    1/(2*sqrt(pi))

    >>> simplify(Ynm(1, -1, theta, phi).expand(func=True))
    sqrt(6)*sqrt(sin(theta)**2)*exp(-I*phi)/(4*sqrt(pi))

    >>> simplify(Ynm(1, 0, theta, phi).expand(func=True))
    sqrt(3)*cos(theta)/(2*sqrt(pi))

    >>> simplify(Ynm(1, 1, theta, phi).expand(func=True))
    -sqrt(6)*sqrt(sin(theta)**2)*exp(I*phi)/(4*sqrt(pi))

    >>> simplify(Ynm(2, -2, theta, phi).expand(func=True))
    sqrt(30)*exp(-2*I*phi)*sin(theta)**2/(8*sqrt(pi))

    >>> simplify(Ynm(2, -1, theta, phi).expand(func=True))
    sqrt(30)*sqrt(sin(theta)**2)*exp(-I*phi)*cos(theta)/(4*sqrt(pi))

    >>> simplify(Ynm(2, 0, theta, phi).expand(func=True))
    sqrt(5)*(-3*sin(theta)**2 + 2)/(4*sqrt(pi))

    >>> simplify(Ynm(2, 1, theta, phi).expand(func=True))
    -sqrt(30)*sqrt(sin(theta)**2)*exp(I*phi)*cos(theta)/(4*sqrt(pi))

    >>> simplify(Ynm(2, 2, theta, phi).expand(func=True))
    sqrt(30)*exp(2*I*phi)*sin(theta)**2/(8*sqrt(pi))

    We can differentiate the functions with respect
    to both angles

    >>> from sympy import Ynm, Symbol, diff
    >>> from sympy.abc import n,m
    >>> theta = Symbol("theta")
    >>> phi = Symbol("phi")

    >>> diff(Ynm(n, m, theta, phi), theta)
    m*cot(theta)*Ynm(n, m, theta, phi) + sqrt((-m + n)*(m + n + 1))*exp(-I*phi)*Ynm(n, m + 1, theta, phi)

    >>> diff(Ynm(n, m, theta, phi), phi)
    I*m*Ynm(n, m, theta, phi)

    Further we can compute the complex conjugation

    >>> from sympy import Ynm, Symbol, conjugate
    >>> from sympy.abc import n,m
    >>> theta = Symbol("theta")
    >>> phi = Symbol("phi")

    >>> conjugate(Ynm(n, m, theta, phi))
    (-1)**(2*m)*exp(-2*I*m*phi)*Ynm(n, m, theta, phi)

    To get back the well known expressions in spherical
    coordinates we use full expansion

    >>> from sympy import Ynm, Symbol, expand_func
    >>> from sympy.abc import n,m
    >>> theta = Symbol("theta")
    >>> phi = Symbol("phi")

    >>> expand_func(Ynm(n, m, theta, phi))
    sqrt((2*n + 1)*(-m + n)!/(m + n)!)*exp(I*m*phi)*assoc_legendre(n, m, cos(theta))/(2*sqrt(pi))

    See Also
    ========

    Ynm_c, Znm

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Spherical_harmonics
    .. [2] http://mathworld.wolfram.com/SphericalHarmonic.html
    .. [3] http://functions.wolfram.com/Polynomials/SphericalHarmonicY/
    """

    nargs = 4

    @classmethod
    def eval(cls, n, m, theta, phi):
        n, m, theta, phi = [sympify(x) for x in (n, m, theta, phi)]

        # Handle negative index m and arguments theta, phi
        if m.could_extract_minus_sign():
            m = -m
            return S.NegativeOne**m * C.exp(-2*I*m*phi) * Ynm(n, m, theta, phi)
        if theta.could_extract_minus_sign():
            theta = -theta
            return Ynm(n, m, theta, phi)
        if phi.could_extract_minus_sign():
            phi = -phi
            return C.exp(-2*I*m*phi) * Ynm(n, m, theta, phi)

        # Move this into the expand_func routine?
        #if n.is_number and m.is_number:
        #    return (sqrt((2*n+1)/(4*pi) * C.factorial(n-m)/C.factorial(n+m)) *
        #            C.exp(I*m*phi) * assoc_legendre(n, m, C.cos(theta)))

        # TODO Add more simplififcation here


    def _eval_expand_func(self, **hints):
        n, m, theta, phi = self.args
        return (sqrt((2*n+1)/(4*pi) * C.factorial(n-m)/C.factorial(n+m)) *
                C.exp(I*m*phi) * assoc_legendre(n, m, C.cos(theta)))

    def fdiff(self, argindex=4):
        if argindex == 1:
            # Diff wrt n
            raise ArgumentIndexError(self, argindex)
        elif argindex == 2:
            # Diff wrt m
            raise ArgumentIndexError(self, argindex)
        elif argindex == 3:
            # Diff wrt theta
            n, m, theta, phi = self.args
            return (m * C.cot(theta) * Ynm(n, m, theta, phi) +
                    sqrt((n-m)*(n+m+1)) * C.exp(-I*phi) * Ynm(n, m+1, theta, phi))
        elif argindex == 4:
            # Diff wrt phi
            n, m, theta, phi = self.args
            return I * m * Ynm(n, m, theta, phi)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_polynomial(self, n, m, theta, phi):
        # TODO: Make sure n \in N
        # TODO: Assert |m| <= n ortherwise we should return 0
        return self.expand(func=True)

    def _eval_conjugate(self):
        # TODO: Make sure theta \in R and phi \in R
        n, m, theta, phi = self.args
        return S.NegativeOne**m * self.func(n, -m, theta, phi)

    def as_real_imag(self, deep=True, **hints):
        # TODO: Handle deep and hints
        n, m, theta, phi = self.args
        re = (sqrt((2*n+1)/(4*pi) * C.factorial(n-m)/C.factorial(n+m)) *
              C.cos(m*phi) * assoc_legendre(n, m, C.cos(theta)))
        im = (sqrt((2*n+1)/(4*pi) * C.factorial(n-m)/C.factorial(n+m)) *
              C.sin(m*phi) * assoc_legendre(n, m, C.cos(theta)))
        return (re, im)


def Ynm_c(n, m, theta, phi):
    r"""Conjugate spherical harmonics defined as

    .. math::
        \overline{Y_n^m(\theta, \varphi)} := (-1)^m Y_n^{-m}(\theta, \varphi)

    See Also
    ========

    Ynm, Znm

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Spherical_harmonics
    .. [2] http://mathworld.wolfram.com/SphericalHarmonic.html
    .. [3] http://functions.wolfram.com/Polynomials/SphericalHarmonicY/
    """
    from sympy import conjugate
    return conjugate(Ynm(n, m, theta, phi))


class Znm(Function):
    r"""
    Real spherical harmonics defined as

    .. math::

        Z_n^m(\theta, \varphi) :=
        \begin{cases}
          \frac{Y_n^m(\theta, \varphi) + \overline{Y_n^m(\theta, \varphi)}}{\sqrt{2}} &\quad m > 0 \\
          Y_n^m(\theta, \varphi) &\quad m = 0 \\
          \frac{Y_n^m(\theta, \varphi) - \overline{Y_n^m(\theta, \varphi)}}{i \sqrt{2}} &\quad m < 0 \\
        \end{cases}

    which gives in simplified form

    .. math::

        Z_n^m(\theta, \varphi) =
        \begin{cases}
          \frac{Y_n^m(\theta, \varphi) + (-1)^m Y_n^{-m}(\theta, \varphi)}{\sqrt{2}} &\quad m > 0 \\
          Y_n^m(\theta, \varphi) &\quad m = 0 \\
          \frac{Y_n^m(\theta, \varphi) - (-1)^m Y_n^{-m}(\theta, \varphi)}{i \sqrt{2}} &\quad m < 0 \\
        \end{cases}

    See Also
    ========

    Ynm, Ynm_c

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Spherical_harmonics
    .. [2] http://mathworld.wolfram.com/SphericalHarmonic.html
    .. [3] http://functions.wolfram.com/Polynomials/SphericalHarmonicY/
    """

    nargs = 4

    @classmethod
    def eval(cls, n, m, theta, phi):
        n, m, th, ph = [sympify(x) for x in (n, m, theta, phi)]

        if m.is_positive:
            zz = (Ynm(n, m, th, ph) + Ynm_c(n, m, th, ph)) / sqrt(2)
            #zz = zz.expand(complex=True)
            #zz = simplify(zz)
            return zz
        elif m.is_zero:
            return Ynm(n, m, th, ph)
        elif m.is_negative:
            zz = (Ynm(n, m, th, ph) - Ynm_c(n, m, th, ph)) / (sqrt(2)*I)
            #zz = zz.expand(complex=True)
            #zz = simplify(zz)
            return zz


def Plmcos(l, m, th):
    """
    Plm(cos(th)).
    """
    l = sympify(l)
    m = sympify(m)
    sin = C.sin
    cos = C.cos
    P = Plm(l, m, _x).subs(_x, cos(th))
    # assume th in (0,pi) => sin(th) is nonnegative
    _sinth = Dummy("_sinth", nonnegative=True)
    P = P.subs(1 - cos(th)**2, _sinth**2).subs(_sinth, sin(th))
    return P


def Ylm(l, m, theta, phi):
    """
    Spherical harmonics Ylm.

    Examples
    ========

    >>> from sympy import symbols, Ylm
    >>> theta, phi = symbols("theta phi")
    >>> Ylm(0, 0, theta, phi)
    1/(2*sqrt(pi))
    >>> Ylm(1, -1, theta, phi)
    sqrt(6)*exp(-I*phi)*sin(theta)/(4*sqrt(pi))
    >>> Ylm(1, 0, theta, phi)
    sqrt(3)*cos(theta)/(2*sqrt(pi))

    """
    l, m, theta, phi = [sympify(x) for x in (l, m, theta, phi)]
    factorial = C.factorial
    return sqrt((2*l + 1)/(4*pi) * factorial(l - m)/factorial(l + m)) * \
        Plmcos(l, m, theta) * C.exp(I*m*phi)


def Ylm_c(l, m, theta, phi):
    """Conjugate spherical harmonics."""
    return (-1)**m * Ylm(l, -m, theta, phi)


def Zlm(l, m, th, ph):
    """
    Real spherical harmonics.
    """
    from sympy import simplify
    if m > 0:
        zz = C.NegativeOne()**m * \
            (Ylm(l, m, th, ph) + Ylm_c(l, m, th, ph))/sqrt(2)
    elif m == 0:
        return Ylm(l, m, th, ph)
    else:
        zz = C.NegativeOne()**m * \
            (Ylm(l, -m, th, ph) - Ylm_c(l, -m, th, ph))/(I*sqrt(2))

    zz = zz.expand(complex=True)
    zz = simplify(zz)
    return zz
