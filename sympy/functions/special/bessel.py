"""Bessel type functions"""

from sympy import S, pi, I
from sympy.core.function import Function, ArgumentIndexError
from sympy.functions.elementary.trigonometric import sin, cos
from sympy.functions.elementary.miscellaneous import sqrt

# TODO
# o Airy Ai and Bi functions
# o Scorer functions G1 and G2
# o Asymptotic expansions
#   These are possible, e.g. for fixed order, but since the bessel type
#   functions are oscillatory they are not actually tractable at
#   infinity, so this is not particularly useful right now.
# o Series Expansions for functions of the second kind about zero
# o Nicer series expansions.
# o More rewriting.
# o Add solvers to ode.py (or rather add solvers for the hypergeometric equation).

class BesselBase(Function):
    """
    Abstract base class for bessel-type functions.

    This class is meant to reduce code duplication.
    All bessel type functions can 1) be differentiated, and the derivatives
    expressed in terms of similar functions and 2) be rewritten in terms
    of other bessel-type functions.

    Here "bessel-type functions" are assumed to have one complex parameter.

    To use this base class, define class attributes ``_a`` and ``_b`` such that
    ``2*F_n' = -_a*F_{n+1} b*F_{n-1}``.
    """

    nargs = 2

    @property
    def order(self):
        """ The order of the bessel-type function. """
        return self.args[0]

    @property
    def argument(self):
        """ The argument of the bessel-type function. """
        return self.args[1]

    def fdiff(self, argindex=2):
        if argindex != 2:
            raise ArgumentIndexError(self, argindex)
        return self._b/2 * self.__class__(self.order - 1, self.argument) \
             - self._a/2 * self.__class__(self.order + 1, self.argument) \

class besselj(BesselBase):
    r"""
    Bessel function of the first kind.

    The Bessel J function of order :math:`\nu` is defined to be the function
    satisfying Bessel's differential equation

    .. math ::
        z^2 \frac{\mathrm{d}^2 w}{\mathrm{d}z^2}
        + z \frac{\mathrm{d}w}{\mathrm{d}z} + (z^2 - \nu^2) w = 0,

    with Laurent expansion

    .. math ::
        J_\nu(z) = z^\nu \left(\frac{1}{\Gamma(\nu + 1) 2^\nu} + O(z^2) \right),

    if :math:`\nu` is not a negative integer. If :math:`\nu=-n \in \mathbb{Z}_{<0}`
    *is* a negative integer, then the definition is

    .. math ::
        J_{-n}(z) = (-1)^n J_n(z).

    Examples
    ========

    Create a bessel function object:

    >>> from sympy import besselj, jn
    >>> from sympy.abc import z, n
    >>> b = besselj(n, z)

    Differentiate it:

    >>> b.diff(z)
    besselj(n - 1, z)/2 - besselj(n + 1, z)/2

    Rewrite in terms of spherical bessel functions:

    >>> b.rewrite(jn)
    sqrt(2)*sqrt(z)*jn(n - 1/2, z)/sqrt(pi)

    Access the parameter and argument:

    >>> b.order
    n
    >>> b.argument
    z

    See Also
    ========

    bessely, besseli, besselk


    References
    ==========

    - Abramowitz, Milton; Stegun, Irene A., eds. (1965), "Chapter 9",
      Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical
      Tables
    - Luke, Y. L. (1969), The Special Functions and Their Approximations,
      Volume 1
    - http://en.wikipedia.org/wiki/Bessel_function
    """

    _a = S.One
    _b = S.One

    def _eval_rewrite_as_jn(self, nu, z, expand=False):
        jn_part = jn(nu - S('1/2'), self.argument)
        if expand:
            jn_part = jn_part._eval_expand_func()
        return sqrt(2*z/pi) * jn_part

    @classmethod
    def eval(cls, nu, z):
        if nu.is_Integer:
            if nu < 0:
                return S(-1)**nu*besselj(-nu, z)
            if z.could_extract_minus_sign():
                return S(-1)**nu*besselj(nu, -z)
            newz = z.extract_multiplicatively(I)
            if newz: # NOTE we don't want to change the function if z==0
                return I**(nu)*besseli(nu, newz)

        # branch handling:
        from sympy import unpolarify, exp
        if nu.is_integer:
            newz = unpolarify(z)
            if newz != z:
                return besselj(nu, newz)
        else:
            newz, n = z.extract_branch_factor()
            if n != 0:
                return exp(2*n*pi*nu*I)*besselj(nu, newz)
        nnu = unpolarify(nu)
        if nu != nnu:
            return besselj(nnu, z)

    def _eval_expand_func(self, **hints):
        if self.order.is_Rational and self.order.q == 2:
            return self._eval_rewrite_as_jn(*self.args, **{'expand':True})
        return self

    def _eval_rewrite_as_besseli(self, nu, z):
        from sympy import polar_lift, exp
        return exp(I*pi*nu/2)*besseli(nu, polar_lift(-I)*z)

class bessely(BesselBase):
    r"""
    Bessel function of the second kind.

    The Bessel Y function of order :math:`\nu` is defined as

    .. math ::
        Y_\nu(z) = \lim_{\mu \to \nu} \frac{J_\mu(z) \cos(\pi \mu)
                                            - J_{-\mu}(z)}{\sin(\pi \mu)},

    where :math:`J_\mu(z)` is the Bessel function of the first kind.

    It is a solution to Bessel's equation, and linearly independent from
    :math:`J_\nu`.

    Examples
    ========

    >>> from sympy import bessely, yn
    >>> from sympy.abc import z, n
    >>> b = bessely(n, z)
    >>> b.diff(z)
    bessely(n - 1, z)/2 - bessely(n + 1, z)/2
    >>> b.rewrite(yn)
    sqrt(2)*sqrt(z)*yn(n - 1/2, z)/sqrt(pi)

    See Also
    ========

    besselj, besseli, besselk

    """

    _a = S.One
    _b = S.One

    def _eval_rewrite_as_yn(self, nu, z, expand=False):
        yn_part = yn(nu - S('1/2'), self.argument)
        if expand:
            yn_part = yn_part._eval_expand_func()
        return sqrt(2*z/pi) * yn_part

    @classmethod
    def eval(cls, nu, z):
        if nu.is_Integer:
            if nu < 0:
                return S(-1)**nu*bessely(-nu, z)

    def _eval_expand_func(self, **hints):
        if self.order.is_Rational and self.order.q == 2:
            return self._eval_rewrite_as_yn(*self.args, **{'expand':True})
        return self

class besseli(BesselBase):
    r"""
    Modified Bessel function of the first kind.

    The Bessel I function is a solution to the modified Bessel equation

    .. math ::
        z^2 \frac{\mathrm{d}^2 w}{\mathrm{d}z^2}
        + z \frac{\mathrm{d}w}{\mathrm{d}z} + (z^2 + \nu^2)^2 w = 0.

    It can be defined as

    .. math ::
        I_\nu(z) = i^{-\nu} J_\nu(iz),

    where :math:`J_\mu(z)` is the Bessel function of the first kind.

    Examples
    ========

    >>> from sympy import besseli
    >>> from sympy.abc import z, n
    >>> besseli(n, z).diff(z)
    besseli(n - 1, z)/2 + besseli(n + 1, z)/2

    See Also
    ========

    besselj, bessely, besselk

    """

    _a = -S.One
    _b = S.One

    @classmethod
    def eval(cls, nu, z):
        if nu.is_Integer:
            newz = z.extract_multiplicatively(I)
            if newz: # NOTE we don't want to change the function if z==0
                return I**(-nu)*besselj(nu, -newz)

        # branch handling:
        from sympy import unpolarify, exp
        if nu.is_integer:
            newz = unpolarify(z)
            if newz != z:
                return besseli(nu, newz)
        else:
            newz, n = z.extract_branch_factor()
            if n != 0:
                return exp(2*n*pi*nu*I)*besseli(nu, newz)
        nnu = unpolarify(nu)
        if nu != nnu:
            return besseli(nnu, z)

    def _eval_rewrite_as_besselj(self, nu, z):
        from sympy import polar_lift, exp
        return exp(-I*pi*nu/2)*besselj(nu, polar_lift(I)*z)

class besselk(BesselBase):
    r"""
    Modified Bessel function of the second kind.

    The Bessel K function of order :math:`\nu` is defined as

    .. math ::
        K_\nu(z) = \lim_{\mu \to \nu} \frac{\pi}{2}
                   \frac{I_{-\mu}(z) -I_\mu(z)}{\sin(\pi \mu)},

    where :math:`I_\mu(z)` is the modified Bessel function of the first kind.

    It is a solution of the modified Bessel equation, and linearly independent
    from :math:`Y_\nu`.

    Examples
    ========

    >>> from sympy import besselk
    >>> from sympy.abc import z, n
    >>> besselk(n, z).diff(z)
    -besselk(n - 1, z)/2 - besselk(n + 1, z)/2

    See Also
    ========

    besselj, besseli, bessely

    """

    _a = S.One
    _b = -S.One

class hankel1(BesselBase):
    r"""
    Hankel function of the first kind.

    This function is defined as

    .. math ::
        H_\nu^{(1)} = J_\nu(z) + iY_\nu(z),

    where :math:`J_\nu(z)` is the Bessel function of the first kind, and
    :math:`Y_\nu(z)` is the Bessel function of the second kind.

    It is a solution to Bessel's equation.

    Examples
    ========

    >>> from sympy import hankel1
    >>> from sympy.abc import z, n
    >>> hankel1(n, z).diff(z)
    hankel1(n - 1, z)/2 - hankel1(n + 1, z)/2

    See Also
    ========

    hankel2, besselj, bessely

    """

    _a = S.One
    _b = S.One

class hankel2(BesselBase):
    r"""
    Hankel function of the second kind.

    This function is defined as

    .. math ::
        H_\nu^{(2)} = J_\nu(z) - iY_\nu(z),

    where :math:`J_\nu(z)` is the Bessel function of the first kind, and
    :math:`Y_\nu(z)` is the Bessel function of the second kind.

    It is a solution to Bessel's equation, and linearly independent from
    :math:`H_\nu^{(1)}`.

    Examples
    ========

    >>> from sympy import hankel2
    >>> from sympy.abc import z, n
    >>> hankel2(n, z).diff(z)
    hankel2(n - 1, z)/2 - hankel2(n + 1, z)/2

    See Also
    ========

    hankel1, besselj, bessely

    """

    _a = S.One
    _b = S.One

from sympy.polys.orthopolys import spherical_bessel_fn as fn

class SphericalBesselBase(BesselBase):
    """
    Base class for spherical bessel functions.

    These are thin wrappers around ordinary bessel functions,
    since spherical bessel functions differ from the ordinary
    ones just by a slight change in order.

    To use this class, define the _rewrite and _expand methods.
    """

    def _expand(self, **hints):
        """ Expand self into a polynomial. Nu is guaranteed to be Integer. """
        raise NotImplementedError('expansion')

    def _rewrite(self):
        """ Rewrite self in terms of ordinary bessel functions. """
        raise NotImplementedError('rewriting')

    def _eval_expand_func(self, **hints):
        if self.order.is_Integer:
            return self._expand(**hints)
        else:
            return self

    def _eval_evalf(self, prec):
        return self._rewrite()._eval_evalf(prec)

    def fdiff(self, argindex=2):
        if argindex != 2:
            raise ArgumentIndexError(self, argindex)
        return self.__class__(self.order - 1, self.argument) - \
               self * (self.order + 1)/self.argument


class jn(SphericalBesselBase):
    r"""
    Spherical Bessel function of the first kind.

    This function is a solution to the spherical bessel equation

    .. math ::
        z^2 \frac{\mathrm{d}^2 w}{\mathrm{d}z^2}
          + 2z \frac{\mathrm{d}w}{\mathrm{d}z} + (z^2 - \nu(\nu + 1)) w = 0.

    It can be defined as

    .. math ::
        j_\nu(z) = \sqrt{\frac{\pi}{2z}} J_{\nu + \frac{1}{2}}(z),

    where :math:`J_\nu(z)` is the Bessel function of the first kind.

    Examples
    ========

    >>> from sympy import Symbol, jn, sin, cos, expand_func
    >>> z = Symbol("z")
    >>> print jn(0, z).expand(func=True)
    sin(z)/z
    >>> jn(1, z).expand(func=True) == sin(z)/z**2 - cos(z)/z
    True
    >>> expand_func(jn(3, z))
    (-6/z**2 + 15/z**4)*sin(z) + (1/z - 15/z**3)*cos(z)

    The spherical Bessel functions of integral order
    are calculated using the formula:

    .. math:: j_n(z) = f_n(z) \sin{z} + (-1)^{n+1} f_{-n-1}(z) \cos{z},

    where the coefficients :math:`f_n(z)` are available as
    :func:`polys.orthopolys.spherical_bessel_fn`.

    See Also
    ========

    besselj, bessely, besselk, yn

    """

    def _rewrite(self):
        return self._eval_rewrite_as_besselj(self.order, self.argument)

    def _eval_rewrite_as_besselj(self, nu, z):
        return sqrt(pi/(2*z)) * besselj(nu + S('1/2'), z)

    def _expand(self, **hints):
        n = self.order
        z = self.argument
        return fn(n, z) * sin(z) + (-1)**(n+1) * fn(-n-1, z) * cos(z)

class yn(SphericalBesselBase):
    r"""
    Spherical Bessel function of the second kind.

    This function is another solution to the spherical bessel equation, and
    linearly independent from :math:`j_n`. It can be defined as

    .. math ::
        j_\nu(z) = \sqrt{\frac{\pi}{2z}} Y_{\nu + \frac{1}{2}}(z),

    where :math:`Y_\nu(z)` is the Bessel function of the second kind.

    Examples
    ========

    >>> from sympy import Symbol, yn, sin, cos, expand_func
    >>> z = Symbol("z")
    >>> print expand_func(yn(0, z))
    -cos(z)/z
    >>> expand_func(yn(1, z)) == -cos(z)/z**2-sin(z)/z
    True

    For integral orders :math:`n`, :math:`y_n` is calculated using the formula:

    .. math:: y_n(z) = (-1)^{n+1} j_{-n-1}(z)

    See Also
    ========

    besselj, bessely, besselk, jn

    """

    def _rewrite(self):
        return self._eval_rewrite_as_bessely(self.order, self.argument)

    def _eval_rewrite_as_bessely(self, nu, z):
        return sqrt(pi/(2*z)) * bessely(nu + S('1/2'), z)

    def _expand(self, **hints):
        n = self.order
        z = self.argument
        return (-1)**(n+1) * \
               (fn(-n-1, z) * sin(z) + (-1)**(-n) * fn(n, z) * cos(z))


def jn_zeros(n, k, method="sympy", dps=15):
    """
    Zeros of the spherical Bessel function of the first kind.

    This returns an array of zeros of jn up to the k-th zero.

    * method = "sympy": uses mpmath besseljzero
    * method = "scipy": uses the SciPy's sph_jn and newton to find all
      roots, which is faster than computing the zeros using a general
      numerical solver, but it requires SciPy and only works with low
      precision floating point numbers.  [the function used with
      method="sympy" is a recent addition to mpmath, before that a general
      solver was used]

    Examples
    ========

    >>> from sympy import jn_zeros
    >>> jn_zeros(2, 4, dps=5)
    [5.7635, 9.095, 12.323, 15.515]

    See Also
    ========

    jn, yn, besselj, besselk, bessely
    """
    from math import pi

    if method == "sympy":
        from sympy.mpmath import besseljzero
        from sympy.mpmath.libmp.libmpf import dps_to_prec
        from sympy import Expr
        prec = dps_to_prec(dps)
        return [Expr._from_mpmath(besseljzero(S(n + 0.5)._to_mpmath(prec),
                                              int(k)), prec) \
                for k in xrange(1, k + 1)]
    elif method == "scipy":
        from scipy.special import sph_jn
        from scipy.optimize import newton
        f  = lambda x: sph_jn(n, x)[0][-1]
    else:
        raise NotImplementedError("Unknown method.")
    def solver(f, x):
        if method == "scipy":
            root = newton(f, x)
        else:
            raise NotImplementedError("Unknown method.")
        return root

    # we need to approximate the position of the first root:
    root = n+pi
    # determine the first root exactly:
    root = solver(f, root)
    roots = [root]
    for i in range(k-1):
        # estimate the position of the next root using the last root + pi:
        root = solver(f, root+pi)
        roots.append(root)
    return roots
