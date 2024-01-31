""" Elliptic Integrals. """
import functools

from sympy.core import S
from sympy.core.function import Function, ArgumentIndexError
from sympy.core.numbers import Rational, Integer, pi, I
from sympy.core.symbol import Dummy
from sympy.functions.elementary.complexes import sign
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.hyperbolic import atanh, sech, tanh
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import cos, sin, tan
from sympy.functions.special.gamma_functions import gamma
from sympy.functions.special.hyper import hyper, meijerg


class elliptic_k(Function):
    r"""
    The complete elliptic integral of the first kind, defined by

    .. math:: K(m) = F\left(\tfrac{\pi}{2}\middle| m\right) =
              \int_0^{\frac{\pi}{2} \frac{dt}{\sqrt{1 - m \sin^2 t}}

    where $F\left(z\middle| m\right)$ is the Legendre incomplete
    elliptic integral of the first kind.

    Explanation
    ===========

    The function $K(m)$ is a single-valued function on the complex
    plane with branch cut along the interval $(1, \infty)$.

    Note that our notation defines the incomplete elliptic integral
    in terms of the parameter $m$ instead of the elliptic modulus
    (eccentricity) $k$.
    In this case, the parameter $m$ is defined as $m=k^2$.

    Examples
    ========

    >>> from sympy import elliptic_k, I
    >>> from sympy.abc import m
    >>> elliptic_k(0)
    pi/2
    >>> elliptic_k(1.0 + I)
    1.50923695405127 + 0.625146415202697*I
    >>> elliptic_k(m).series(n=3)
    pi/2 + pi*m/8 + 9*pi*m**2/128 + O(m**3)

    See Also
    ========

    elliptic_f

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Elliptic_integrals
    .. [2] https://functions.wolfram.com/EllipticIntegrals/EllipticK

    """

    @classmethod
    def eval(cls, m):
        if m.is_zero:
            return pi*S.Half
        elif m is S.Half:
            return 8*pi**Rational(3, 2)/gamma(Rational(-1, 4))**2
        elif m is S.One:
            return S.ComplexInfinity
        elif m is S.NegativeOne:
            return gamma(Rational(1, 4))**2/(4*sqrt(2*pi))
        elif m in (S.Infinity, S.NegativeInfinity, I*S.Infinity,
                   I*S.NegativeInfinity, S.ComplexInfinity):
            return S.Zero

    def fdiff(self, argindex=1):
        m = self.args[0]
        return (elliptic_e(m) - (1 - m)*elliptic_k(m))/(2*m*(1 - m))

    def _eval_conjugate(self):
        m = self.args[0]
        if (m.is_real and (m - 1).is_positive) is False:
            return self.func(m.conjugate())

    def _eval_nseries(self, x, n, logx, cdir=0):
        from sympy.simplify import hyperexpand
        return hyperexpand(self.rewrite(hyper)._eval_nseries(x, n=n, logx=logx))

    def _eval_rewrite_as_hyper(self, m, **kwargs):
        return pi*S.Half*hyper((S.Half, S.Half), (S.One,), m)

    def _eval_rewrite_as_meijerg(self, m, **kwargs):
        return meijerg(((S.Half, S.Half), []), ((S.Zero,), (S.Zero,)), -m)/2

    def _eval_is_zero(self):
        m = self.args[0]
        if m.is_infinite:
            return True

    def _eval_rewrite_as_Integral(self, *args, **kwargs):
        from sympy.integrals.integrals import Integral
        t = Dummy('t')
        m = self.args[0]
        return Integral(1/sqrt(1 - m*sin(t)**2), (t, 0, pi/2))


class elliptic_f(Function):
    r"""
    The Legendre incomplete elliptic integral of the first
    kind, defined by

    .. math:: F\left(z\middle| m\right) =
              \int_0^z \frac{dt}{\sqrt{1 - m \sin^2 t}}

    Explanation
    ===========

    This function reduces to a complete elliptic integral of
    the first kind, $K(m)$, when $z = \pi/2$.

    Note that our notation defines the incomplete elliptic integral
    in terms of the parameter $m$ instead of the elliptic modulus
    (eccentricity) $k$.
    In this case, the parameter $m$ is defined as $m=k^2$.

    Examples
    ========

    >>> from sympy import elliptic_f, I
    >>> from sympy.abc import z, m
    >>> elliptic_f(z, m).series(z)
    z + z**5*(3*m**2/40 - m/30) + m*z**3/6 + O(z**6)
    >>> elliptic_f(3.0 + I/2, 1.0 + I)
    2.909449841483 + 1.74720545502474*I

    See Also
    ========

    elliptic_k

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Elliptic_integrals
    .. [2] https://functions.wolfram.com/EllipticIntegrals/EllipticF

    """

    @classmethod
    def eval(cls, z, m):
        if z.is_zero:
            return S.Zero
        if m.is_zero:
            return z
        k = 2*z/pi
        if k.is_integer:
            return k*elliptic_k(m)
        elif m in (S.Infinity, S.NegativeInfinity):
            return S.Zero
        elif z.could_extract_minus_sign():
            return -elliptic_f(-z, m)

    def fdiff(self, argindex=1):
        z, m = self.args
        fm = sqrt(1 - m*sin(z)**2)
        if argindex == 1:
            return 1/fm
        elif argindex == 2:
            return (elliptic_e(z, m)/(2*m*(1 - m)) - elliptic_f(z, m)/(2*m) -
                    sin(2*z)/(4*(1 - m)*fm))
        raise ArgumentIndexError(self, argindex)

    def _eval_conjugate(self):
        z, m = self.args
        if (m.is_real and (m - 1).is_positive) is False:
            return self.func(z.conjugate(), m.conjugate())

    def _eval_rewrite_as_Integral(self, *args, **kwargs):
        from sympy.integrals.integrals import Integral
        t = Dummy('t')
        z, m = self.args[0], self.args[1]
        return Integral(1/(sqrt(1 - m*sin(t)**2)), (t, 0, z))

    def _eval_is_zero(self):
        z, m = self.args
        if z.is_zero:
            return True
        if m.is_extended_real and m.is_infinite:
            return True


class elliptic_e(Function):
    r"""
    Called with two arguments $z$ and $m$, evaluates the
    incomplete elliptic integral of the second kind, defined by

    .. math:: E\left(z\middle| m\right) = \int_0^z \sqrt{1 - m \sin^2 t} dt

    Called with a single argument $m$, evaluates the Legendre complete
    elliptic integral of the second kind

    .. math:: E(m) = E\left(\tfrac{\pi}{2}\middle| m\right)

    Explanation
    ===========

    The function $E(m)$ is a single-valued function on the complex
    plane with branch cut along the interval $(1, \infty)$.

    Note that our notation defines the incomplete elliptic integral
    in terms of the parameter $m$ instead of the elliptic modulus
    (eccentricity) $k$.
    In this case, the parameter $m$ is defined as $m=k^2$.

    Examples
    ========

    >>> from sympy import elliptic_e, I
    >>> from sympy.abc import z, m
    >>> elliptic_e(z, m).series(z)
    z + z**5*(-m**2/40 + m/30) - m*z**3/6 + O(z**6)
    >>> elliptic_e(m).series(n=4)
    pi/2 - pi*m/8 - 3*pi*m**2/128 - 5*pi*m**3/512 + O(m**4)
    >>> elliptic_e(1 + I, 2 - I/2).n()
    1.55203744279187 + 0.290764986058437*I
    >>> elliptic_e(0)
    pi/2
    >>> elliptic_e(2.0 - I)
    0.991052601328069 + 0.81879421395609*I

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Elliptic_integrals
    .. [2] https://functions.wolfram.com/EllipticIntegrals/EllipticE2
    .. [3] https://functions.wolfram.com/EllipticIntegrals/EllipticE

    """

    @classmethod
    def eval(cls, m, z=None):
        if z is not None:
            z, m = m, z
            k = 2*z/pi
            if m.is_zero:
                return z
            if z.is_zero:
                return S.Zero
            elif k.is_integer:
                return k*elliptic_e(m)
            elif m in (S.Infinity, S.NegativeInfinity):
                return S.ComplexInfinity
            elif z.could_extract_minus_sign():
                return -elliptic_e(-z, m)
        else:
            if m.is_zero:
                return pi/2
            elif m is S.One:
                return S.One
            elif m is S.Infinity:
                return I*S.Infinity
            elif m is S.NegativeInfinity:
                return S.Infinity
            elif m is S.ComplexInfinity:
                return S.ComplexInfinity

    def fdiff(self, argindex=1):
        if len(self.args) == 2:
            z, m = self.args
            if argindex == 1:
                return sqrt(1 - m*sin(z)**2)
            elif argindex == 2:
                return (elliptic_e(z, m) - elliptic_f(z, m))/(2*m)
        else:
            m = self.args[0]
            if argindex == 1:
                return (elliptic_e(m) - elliptic_k(m))/(2*m)
        raise ArgumentIndexError(self, argindex)

    def _eval_conjugate(self):
        if len(self.args) == 2:
            z, m = self.args
            if (m.is_real and (m - 1).is_positive) is False:
                return self.func(z.conjugate(), m.conjugate())
        else:
            m = self.args[0]
            if (m.is_real and (m - 1).is_positive) is False:
                return self.func(m.conjugate())

    def _eval_nseries(self, x, n, logx, cdir=0):
        from sympy.simplify import hyperexpand
        if len(self.args) == 1:
            return hyperexpand(self.rewrite(hyper)._eval_nseries(x, n=n, logx=logx))
        return super()._eval_nseries(x, n=n, logx=logx)

    def _eval_rewrite_as_hyper(self, *args, **kwargs):
        if len(args) == 1:
            m = args[0]
            return (pi/2)*hyper((Rational(-1, 2), S.Half), (S.One,), m)

    def _eval_rewrite_as_meijerg(self, *args, **kwargs):
        if len(args) == 1:
            m = args[0]
            return -meijerg(((S.Half, Rational(3, 2)), []), \
                            ((S.Zero,), (S.Zero,)), -m)/4

    def _eval_rewrite_as_Integral(self, *args, **kwargs):
        from sympy.integrals.integrals import Integral
        z, m = (pi/2, self.args[0]) if len(self.args) == 1 else self.args
        t = Dummy('t')
        return Integral(sqrt(1 - m*sin(t)**2), (t, 0, z))


class elliptic_pi(Function):
    r"""
    Called with three arguments $n$, $z$ and $m$, evaluates the
    Legendre incomplete elliptic integral of the third kind, defined by

    .. math:: \Pi\left(n; z\middle| m\right) = \int_0^z \frac{dt}
              {\left(1 - n \sin^2 t\right) \sqrt{1 - m \sin^2 t}}

    Called with two arguments $n$ and $m$, evaluates the complete
    elliptic integral of the third kind:

    .. math:: \Pi\left(n\middle| m\right) =
              \Pi\left(n; \tfrac{\pi}{2}\middle| m\right)

    Explanation
    ===========

    Note that our notation defines the incomplete elliptic integral
    in terms of the parameter $m$ instead of the elliptic modulus
    (eccentricity) $k$.
    In this case, the parameter $m$ is defined as $m=k^2$.

    Examples
    ========

    >>> from sympy import elliptic_pi, I
    >>> from sympy.abc import z, n, m
    >>> elliptic_pi(n, z, m).series(z, n=4)
    z + z**3*(m/6 + n/3) + O(z**4)
    >>> elliptic_pi(0.5 + I, 1.0 - I, 1.2)
    2.50232379629182 - 0.760939574180767*I
    >>> elliptic_pi(0, 0)
    pi/2
    >>> elliptic_pi(1.0 - I/3, 2.0 + I)
    3.29136443417283 + 0.32555634906645*I

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Elliptic_integrals
    .. [2] https://functions.wolfram.com/EllipticIntegrals/EllipticPi3
    .. [3] https://functions.wolfram.com/EllipticIntegrals/EllipticPi

    """

    @classmethod
    def eval(cls, n, m, z=None):
        if z is not None:
            n, z, m = n, m, z
            if n.is_zero:
                return elliptic_f(z, m)
            elif n is S.One:
                return (elliptic_f(z, m) +
                        (sqrt(1 - m*sin(z)**2)*tan(z) -
                         elliptic_e(z, m))/(1 - m))
            k = 2*z/pi
            if k.is_integer:
                return k*elliptic_pi(n, m)
            elif m.is_zero:
                return atanh(sqrt(n - 1)*tan(z))/sqrt(n - 1)
            elif n == m:
                return (elliptic_f(z, n) - elliptic_pi(1, z, n) +
                        tan(z)/sqrt(1 - n*sin(z)**2))
            elif n in (S.Infinity, S.NegativeInfinity):
                return S.Zero
            elif m in (S.Infinity, S.NegativeInfinity):
                return S.Zero
            elif z.could_extract_minus_sign():
                return -elliptic_pi(n, -z, m)
            if n.is_zero:
                return elliptic_f(z, m)
            if m.is_extended_real and m.is_infinite or \
                    n.is_extended_real and n.is_infinite:
                return S.Zero
        else:
            if n.is_zero:
                return elliptic_k(m)
            elif n is S.One:
                return S.ComplexInfinity
            elif m.is_zero:
                return pi/(2*sqrt(1 - n))
            elif m == S.One:
                return S.NegativeInfinity/sign(n - 1)
            elif n == m:
                return elliptic_e(n)/(1 - n)
            elif n in (S.Infinity, S.NegativeInfinity):
                return S.Zero
            elif m in (S.Infinity, S.NegativeInfinity):
                return S.Zero
            if n.is_zero:
                return elliptic_k(m)
            if m.is_extended_real and m.is_infinite or \
                    n.is_extended_real and n.is_infinite:
                return S.Zero

    def _eval_conjugate(self):
        if len(self.args) == 3:
            n, z, m = self.args
            if (n.is_real and (n - 1).is_positive) is False and \
               (m.is_real and (m - 1).is_positive) is False:
                return self.func(n.conjugate(), z.conjugate(), m.conjugate())
        else:
            n, m = self.args
            return self.func(n.conjugate(), m.conjugate())

    def fdiff(self, argindex=1):
        if len(self.args) == 3:
            n, z, m = self.args
            fm, fn = sqrt(1 - m*sin(z)**2), 1 - n*sin(z)**2
            if argindex == 1:
                return (elliptic_e(z, m) + (m - n)*elliptic_f(z, m)/n +
                        (n**2 - m)*elliptic_pi(n, z, m)/n -
                        n*fm*sin(2*z)/(2*fn))/(2*(m - n)*(n - 1))
            elif argindex == 2:
                return 1/(fm*fn)
            elif argindex == 3:
                return (elliptic_e(z, m)/(m - 1) +
                        elliptic_pi(n, z, m) -
                        m*sin(2*z)/(2*(m - 1)*fm))/(2*(n - m))
        else:
            n, m = self.args
            if argindex == 1:
                return (elliptic_e(m) + (m - n)*elliptic_k(m)/n +
                        (n**2 - m)*elliptic_pi(n, m)/n)/(2*(m - n)*(n - 1))
            elif argindex == 2:
                return (elliptic_e(m)/(m - 1) + elliptic_pi(n, m))/(2*(n - m))
        raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_Integral(self, *args, **kwargs):
        from sympy.integrals.integrals import Integral
        if len(self.args) == 2:
            n, m, z = self.args[0], self.args[1], pi/2
        else:
            n, z, m = self.args
        t = Dummy('t')
        return Integral(1/((1 - n*sin(t)**2)*sqrt(1 - m*sin(t)**2)), (t, 0, z))


class ThetaBase(Function):
    r"""
    Base class for Jacobi theta function $\vartheta_n(z,q)$, where $n=1,2,3,4$.

    See Also
    ========

    theta1, theta2, theta3, theta4

    """

    def __new__(cls, z, q=0):
        if q is None:
            return super().__new__(cls, 0, z)

        return super().__new__(cls, z, q)

    def _latex(self, printer):
        args = ", ".join(printer._print(arg) for arg in self.args)
        return r"\vartheta_{}\left({}\right)".format(self._type_val, args)

    def _eval_mpmath(self):
        from mpmath import jtheta
        return jtheta, (Integer(self._type_val), *self.args)


class theta1(ThetaBase):
    r"""
    Jacobi theta function $\vartheta_1(z,q)$ defined as

    .. math:: \vartheta_1(z,q) = 2 q^{\sfrac{1}{4}} \sum_{n=0}^\infty
              (-1)^n q^{n^2 + n} \sin((2n + 1)z)

    See Also
    ========

    theta2, theta3, theta4

    """

    _type_val = 1

    @classmethod
    def eval(cls, z, q):
        if q.is_zero or z.is_zero or (z/pi).is_integer:
            return S.Zero

    def _eval_rewrite_as_theta2(self, *args):
        z, q = self.args
        return theta2(z - pi/2, q)

    def _eval_rewrite_as_theta3(self, *args):
        z, q = self.args
        return -I*exp(I*z)*q**(Rational(1, 4))*theta3(z + (pi - I*log(q))/2, q)

    def _eval_rewrite_as_theta4(self, *args):
        z, q = self.args
        return I*exp(-I*z)*q**(Rational(1, 4))*theta4(z + I*log(q)/2, q)


class theta2(ThetaBase):
    r"""
    Jacobi theta function $\vartheta_2(z,q)$ defined as

    .. math:: \vartheta_2(z,q) = 2 q^{\sfrac{1}{4}} \sum_{n=0}^\infty
              q^{n^2 + n} \cos((2n + 1)z)

    See Also
    ========

    theta1, theta3, theta4
    """

    _type_val = 2

    @classmethod
    def eval(cls, z, q):
        if q.is_zero or z in (pi/2, -pi/2) or ((z - pi/2)/pi).is_integer:
            return S.Zero

    def _eval_rewrite_as_theta1(self, *args):
        z, q = self.args
        return theta1(z + pi/2, q)

    def _eval_rewrite_as_theta3(self, *args):
        z, q = self.args
        return exp(I*z)*q**(Rational(1, 4))*theta3(z - I*log(q)/2, q)

    def _eval_rewrite_as_theta4(self, *args):
        z, q = self.args
        return exp(I*z)*q**(Rational(1, 4))*theta4(z + (pi - I*log(q))/2, q)


class theta3(ThetaBase):
    r"""
    Jacobi theta function $\vartheta_3(z,q)$.

    .. math:: \vartheta_3(z,q) = 1 + 2 \sum_{n=0}^\infty q^{n^2 + n} \cos(2nz)

    See Also
    ========

    theta1, theta2, theta4
    """

    _type_val = 3

    @classmethod
    def eval(cls, z, q):
        if q.is_zero:
            return S.One
        if z.is_zero and q.could_extract_minus_sign():
            return theta4(z, -q)

    def _eval_rewrite_as_theta1(self, *args):
        z, q = self.args
        return exp(-I*z)*q**(Rational(1, 4))*theta1(z + (pi + I*log(q))/2, q)

    def _eval_rewrite_as_theta2(self, *args):
        z, q = self.args
        return exp(I*z)*q**(Rational(1, 4))*theta2(z - I*log(q)/2, q)

    def _eval_rewrite_as_theta4(self, *args):
        z, q = self.args
        return theta4(z + pi/2, q)


class theta4(ThetaBase):
    r"""
    Jacobi theta function $\vartheta_4(z,q)$.

    .. math:: \vartheta_4(z,q) = 1 + 2 \sum_{n=0}^\infty (-q)^{n^2 + n} \cos(2nz)

    See Also
    ========

    theta1, theta2, theta4
    """

    _type_val = 4

    @classmethod
    def eval(cls, z, q):
        if q.is_zero:
            return S.One
        if z.is_zero and q.could_extract_minus_sign():
            return theta3(z, -q)

    def _eval_rewrite_as_theta1(self, *args):
        z, q = self.args
        return I*exp(-I*z)*q**(Rational(1, 4))*theta1(z + I*log(q)/2, q)

    def _eval_rewrite_as_theta2(self, *args):
        z, q = self.args
        return I*exp(-I*z)*q**(Rational(1, 4))*theta2(z - (pi - I*log(q))/2, q)

    def _eval_rewrite_as_theta3(self, *args):
        z, q = self.args
        return theta3(z + pi/2, q)


class JacobiEllipticFunctionBase(Function):
    r"""
    Base class for Jacobi elliptic functions.

    See also
    ========

    jacobi_cd, jacobi_cn, jacobi_cs, jacobi_dc, jacobi_dn, jacobi_ds,
    jacobi_nc, jacobi_nd, jacobi_ns, jacobi_sc, jacobi_sd, jacobi_sn

    """
    def _eval_mpmath(self):
        from mpmath import ellipfun
        ellipfunspec = functools.partial(ellipfun, self._type_str)
        return ellipfunspec, self.args

    def _latex(self, printer):
        args = ", ".join(printer._print(arg) for arg in self.args)
        return r"\operatorname{{{}}}\left({}\right)".format(self._type_str, args)

    @staticmethod
    def _eval(type_str, u, m):
        c1, c2 = str(type_str)
        utmp = Dummy()
        u0vals = {'c': S.One, 'd': S.One, 'n': S.One, 's': S.Zero}
        m0vals = {'c': cos(utmp), 'd': S.One, 'n': S.One, 's': sin(utmp)}
        m1vals = {'c': sech(utmp), 'd': sech(utmp), 'n': S.One, 's': tanh(utmp)}
        if u is S.Zero:
            return u0vals[c1]/u0vals[c2]
        if m is S.Zero:
            return (m0vals[c1]/m0vals[c2]).subs(utmp, u)
        if m is S.One:
            return (m1vals[c1]/m1vals[c2]).subs(utmp, u)

    def doit(self, **hints):
        return self._eval_doit(**hints)

    def _eval_doit(self, **hints):
        u, m = self.args
        if hints.get('deep', True):
            u, m = u.doit(**hints), m.doit(**hints)
        c1, c2 = str(self._type_str)
        utmp = Dummy()
        u0vals = {'c': S.One, 'd': S.One, 'n': S.One, 's': S.Zero}
        m0vals = {'c': cos(utmp), 'd': S.One, 'n': S.One, 's': sin(utmp)}
        m1vals = {'c': sech(utmp), 'd': sech(utmp), 'n': S.One, 's': tanh(utmp)}
        if u.is_zero:
            return u0vals[c1]/u0vals[c2]
        if m.is_zero:
            return (m0vals[c1]/m0vals[c2]).subs(utmp, u)
        if m.equals(S.One):
            return (m1vals[c1]/m1vals[c2]).subs(utmp, u)

    def _eval_expand_func(self, **kwargs):
        num, den = self._type_str
        u, m = self.args
        if den == 'n':
            return self
        funcs = {'c': jacobi_cn, 'd': jacobi_dn, 'n': jacobi_sn}
        if num == 'n':
            return 1/funcs[den](u, m)
        return funcs[num](u, m)/funcs[den](u, m)

    def _eval_rewrite_as_theta2(self, *args):
        expanded = self._eval_expand_func()
        if expanded.has(jacobi_cn, jacobi_sn, jacobi_dn):
            return expanded.rewrite(theta2)
        return self


class jacobi_cd(JacobiEllipticFunctionBase):
    r"""
    Jacobi elliptic function:

    .. math::

        \operatorname{cd}(u \mid m) = \frac{\cos(\phi)}{\sqrt{1 - m \sin^2(\phi)}

    where $\phi = \operatorname{am}(u \mid m)$.

    See also
    ========

    jacobi_cd, jacobi_cn, jacobi_cs, jacobi_dc, jacobi_dn, jacobi_ds,
    jacobi_nc, jacobi_nd, jacobi_ns, jacobi_sc, jacobi_sd, jacobi_sn

    """
    _type_str = "cd"

    @classmethod
    def eval(cls, u, m):
        if u.is_Number and m.is_Number:
            return JacobiEllipticFunctionBase._eval("cd", u, m)

    def fdiff(self, argindex=1):
        u, m = self.args
        if argindex == 1:
            return (m - 1)*jacobi_nd(u, m)*jacobi_sd(u, m)
        if argindex == 2:
            raise NotImplementedError
        raise ArgumentIndexError(self, argindex)


class jacobi_cn(JacobiEllipticFunctionBase):
    r"""
    Jacobi elliptic function: (wrong)

    .. math::

        \operatorname{cd}(u \mid m) = \frac{\cos(\phi)}{\sqrt{1 - m \sin^2(\phi)}

    where $\phi = \operatorname{am}(u \mid m)$.

    See also
    ========

    jacobi_cd, jacobi_cn, jacobi_cs, jacobi_dc, jacobi_dn, jacobi_ds,
    jacobi_nc, jacobi_nd, jacobi_ns, jacobi_sc, jacobi_sd, jacobi_sn

    """
    _type_str = "cn"

    @classmethod
    def eval(cls, u, m):
        if u.is_Number and m.is_Number:
            return JacobiEllipticFunctionBase._eval("cm", u, m)

    def fdiff(self, argindex=1):
        u, m = self.args
        if argindex == 1:
            return -jacobi_sn(u, m)*jacobi_dn(u, m)
        if argindex == 2:
            raise NotImplementedError
        raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_theta2(self, *args):
        u, m = self.args
        q = elliptic_nome_q(m)
        t = u/theta3(0, q)**2
        return theta4(0, q)/theta2(0, q)*theta2(t, q)/theta4(t, q)


class jacobi_cs(JacobiEllipticFunctionBase):
    r"""
    Jacobi elliptic function: (wrong)

    .. math::

        \operatorname{cd}(u \mid m) = \frac{\cos(\phi)}{\sqrt{1 - m \sin^2(\phi)}

    where $\phi = \operatorname{am}(u \mid m)$.

    See also
    ========

    jacobi_cd, jacobi_cn, jacobi_cs, jacobi_dc, jacobi_dn, jacobi_ds,
    jacobi_nc, jacobi_nd, jacobi_ns, jacobi_sc, jacobi_sd, jacobi_sn

    """
    _type_str = "cs"

    @classmethod
    def eval(cls, u, m):
        if u.is_Number and m.is_Number:
            return JacobiEllipticFunctionBase._eval("cs", u, m)

    def fdiff(self, argindex=1):
        u, m = self.args
        if argindex == 1:
            return -jacobi_ds(u, m)*jacobi_ns(u, m)
        if argindex == 2:
            raise NotImplementedError
        raise ArgumentIndexError(self, argindex)


class jacobi_dc(JacobiEllipticFunctionBase):
    r"""
    Jacobi elliptic function:

    .. math::

        \operatorname{cd}(u \mid m) = \frac{\cos(\phi)}{\sqrt{1 - m \sin^2(\phi)}

    where $\phi = \operatorname{am}(u \mid m)$.

    See also
    ========

    jacobi_cd, jacobi_cn, jacobi_cs, jacobi_dc, jacobi_dn, jacobi_ds,
    jacobi_nc, jacobi_nd, jacobi_ns, jacobi_sc, jacobi_sd, jacobi_sn

    """
    _type_str = "dc"

    @classmethod
    def eval(cls, u, m):
        if u.is_Number and m.is_Number:
            return JacobiEllipticFunctionBase._eval("dc", u, m)

    def fdiff(self, argindex=1):
        u, m = self.args
        if argindex == 1:
            return (1 - m)*jacobi_nc(u, m)*jacobi_sd(u, m)
        if argindex == 2:
            raise NotImplementedError
        raise ArgumentIndexError(self, argindex)


class jacobi_dn(JacobiEllipticFunctionBase):
    r"""
    Jacobi elliptic function: (wrong)

    .. math::

        \operatorname{cd}(u \mid m) = \frac{\cos(\phi)}{\sqrt{1 - m \sin^2(\phi)}

    where $\phi = \operatorname{am}(u \mid m)$.

    See also
    ========

    jacobi_cd, jacobi_cn, jacobi_cs, jacobi_dc, jacobi_dn, jacobi_ds,
    jacobi_nc, jacobi_nd, jacobi_ns, jacobi_sc, jacobi_sd, jacobi_sn

    """
    _type_str = "dn"

    @classmethod
    def eval(cls, u, m):
        if u.is_Number and m.is_Number:
            return JacobiEllipticFunctionBase._eval("dm", u, m)

    def fdiff(self, argindex=1):
        u, m = self.args
        if argindex == 1:
            return -m*jacobi_sn(u, m)*jacobi_cn(u, m)
        if argindex == 2:
            raise NotImplementedError
        raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_theta2(self, *args):
        u, m = self.args
        q = elliptic_nome_q(m)
        t = u/theta3(0, q)**2
        return theta4(0, q)/theta3(0, q)*theta3(t, q)/theta4(t, q)


class jacobi_ds(JacobiEllipticFunctionBase):
    r"""
    Jacobi elliptic function: (wrong)

    .. math::

        \operatorname{cd}(u \mid m) = \frac{\cos(\phi)}{\sqrt{1 - m \sin^2(\phi)}

    where $\phi = \operatorname{am}(u \mid m)$.

    See also
    ========

    jacobi_cd, jacobi_cn, jacobi_cs, jacobi_dc, jacobi_dn, jacobi_ds,
    jacobi_nc, jacobi_nd, jacobi_ns, jacobi_sc, jacobi_sd, jacobi_sn

    """
    _type_str = "ds"

    @classmethod
    def eval(cls, u, m):
        if u.is_Number and m.is_Number:
            return JacobiEllipticFunctionBase._eval("ds", u, m)

    def fdiff(self, argindex=1):
        u, m = self.args
        if argindex == 1:
            return -jacobi_cs(u, m)*jacobi_ns(u, m)
        if argindex == 2:
            raise NotImplementedError
        raise ArgumentIndexError(self, argindex)


class jacobi_nc(JacobiEllipticFunctionBase):
    r"""
    Jacobi elliptic function: (wrong)

    .. math::

        \operatorname{cd}(u \mid m) = \frac{\cos(\phi)}{\sqrt{1 - m \sin^2(\phi)}

    where $\phi = \operatorname{am}(u \mid m)$.

    See also
    ========

    jacobi_cd, jacobi_cn, jacobi_cs, jacobi_dc, jacobi_dn, jacobi_ds,
    jacobi_nc, jacobi_nd, jacobi_ns, jacobi_sc, jacobi_sd, jacobi_sn

    """
    _type_str = "nc"

    @classmethod
    def eval(cls, u, m):
        if u.is_Number and m.is_Number:
            return JacobiEllipticFunctionBase._eval("nc", u, m)

    def fdiff(self, argindex=1):
        u, m = self.args
        if argindex == 1:
            return jacobi_dc(u, m)*jacobi_sc(u, m)
        if argindex == 2:
            raise NotImplementedError
        raise ArgumentIndexError(self, argindex)


class jacobi_nd(JacobiEllipticFunctionBase):
    r"""
    Jacobi elliptic function: (wrong)

    .. math::

        \operatorname{cd}(u \mid m) = \frac{\cos(\phi)}{\sqrt{1 - m \sin^2(\phi)}

    where $\phi = \operatorname{am}(u \mid m)$.

    See also
    ========

    jacobi_cd, jacobi_cn, jacobi_cs, jacobi_dc, jacobi_dn, jacobi_ds,
    jacobi_nc, jacobi_nd, jacobi_ns, jacobi_sc, jacobi_sd, jacobi_sn

    """
    _type_str = "nd"

    @classmethod
    def eval(cls, u, m):
        if u.is_Number and m.is_Number:
            return JacobiEllipticFunctionBase._eval("nd", u, m)

    def fdiff(self, argindex=1):
        u, m = self.args
        if argindex == 1:
            return m*jacobi_cd(u, m)*jacobi_sd(u, m)
        if argindex == 2:
            raise NotImplementedError
        raise ArgumentIndexError(self, argindex)


class jacobi_ns(JacobiEllipticFunctionBase):
    r"""
    Jacobi elliptic function: (wrong)

    .. math::

        \operatorname{cd}(u \mid m) = \frac{\cos(\phi)}{\sqrt{1 - m \sin^2(\phi)}

    where $\phi = \operatorname{am}(u \mid m)$.

    See also
    ========

    jacobi_cd, jacobi_cn, jacobi_cs, jacobi_dc, jacobi_dn, jacobi_ds,
    jacobi_nc, jacobi_nd, jacobi_ns, jacobi_sc, jacobi_sd, jacobi_sn

    """
    _type_str = "ns"

    @classmethod
    def eval(cls, u, m):
        if u.is_Number and m.is_Number:
            return JacobiEllipticFunctionBase._eval("ns", u, m)

    def fdiff(self, argindex=1):
        u, m = self.args
        if argindex == 1:
            return -jacobi_cs(u, m)*jacobi_ds(u, m)
        if argindex == 2:
            raise NotImplementedError
        raise ArgumentIndexError(self, argindex)


class jacobi_sc(JacobiEllipticFunctionBase):
    r"""
    Jacobi elliptic function: (wrong)

    .. math::

        \operatorname{cd}(u \mid m) = \frac{\cos(\phi)}{\sqrt{1 - m \sin^2(\phi)}

    where $\phi = \operatorname{am}(u \mid m)$.

    See also
    ========

    jacobi_cd, jacobi_cn, jacobi_cs, jacobi_dc, jacobi_dn, jacobi_ds,
    jacobi_nc, jacobi_nd, jacobi_ns, jacobi_sc, jacobi_sd, jacobi_sn

    """
    _type_str = "sc"

    @classmethod
    def eval(cls, u, m):
        if u.is_Number and m.is_Number:
            return JacobiEllipticFunctionBase._eval("sc", u, m)

    def fdiff(self, argindex=1):
        u, m = self.args
        if argindex == 1:
            return jacobi_dc(u, m)*jacobi_nc(u, m)
        if argindex == 2:
            raise NotImplementedError
        raise ArgumentIndexError(self, argindex)


class jacobi_sd(JacobiEllipticFunctionBase):
    r"""
    Jacobi elliptic function: (wrong)

    .. math::

        \operatorname{cd}(u \mid m) = \frac{\cos(\phi)}{\sqrt{1 - m \sin^2(\phi)}

    where $\phi = \operatorname{am}(u \mid m)$.

    See also
    ========

    jacobi_cd, jacobi_cn, jacobi_cs, jacobi_dc, jacobi_dn, jacobi_ds,
    jacobi_nc, jacobi_nd, jacobi_ns, jacobi_sc, jacobi_sd, jacobi_sn

    """
    _type_str = "sd"

    @classmethod
    def eval(cls, u, m):
        if u.is_Number and m.is_Number:
            return JacobiEllipticFunctionBase._eval("sd", u, m)

    def fdiff(self, argindex=1):
        u, m = self.args
        if argindex == 1:
            return jacobi_cd(u, m)*jacobi_nd(u, m)
        if argindex == 2:
            raise NotImplementedError
        raise ArgumentIndexError(self, argindex)


class jacobi_sn(JacobiEllipticFunctionBase):
    r"""
    Jacobi elliptic function: (wrong)

    .. math::

        \operatorname{cd}(u \mid m) = \frac{\cos(\phi)}{\sqrt{1 - m \sin^2(\phi)}

    where $\phi = \operatorname{am}(u \mid m)$.

    See also
    ========

    jacobi_cd, jacobi_cn, jacobi_cs, jacobi_dc, jacobi_dn, jacobi_ds,
    jacobi_nc, jacobi_nd, jacobi_ns, jacobi_sc, jacobi_sd, jacobi_sn

    """
    _type_str = "sn"

    @classmethod
    def eval(cls, u, m):
        if u.is_Number and m.is_Number:
            return JacobiEllipticFunctionBase._eval("sn", u, m)

    def fdiff(self, argindex=1):
        u, m = self.args
        if argindex == 1:
            return -jacobi_cn(u, m)*jacobi_dn(u, m)
        if argindex == 2:
            raise NotImplementedError
        raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_theta2(self, *args):
        u, m = self.args
        q = elliptic_nome_q(m)
        t = u/theta3(0, q)**2
        return theta3(0, q)/theta2(0, q)*theta3(1, q)/theta4(t, q)


class elliptic_nome_q(Function):
    r"""
    Elliptic nome $q(m)$, defined as

    .. math::

        q(m) = \exp\left(-\frac{\pi K(1-m)}{K(m)}\right)

    where $K(m)$ is :class:`elliptic_k`.

    """
    @classmethod
    def eval(cls, m):
        if m.is_Number:
            if m is S.Zero:
                return S.Zero
            if m is S.One:
                return S.One
            if m is S.Half:
                return exp(-pi)

    def doit(self, **hints):
        return self._eval_doit(**hints)

    def _latex(self, printer):
        return fr"q\left({self.args[0]}\right)"

    def _eval_doit(self, **hints):
        m = self.args[0]
        if hints.get('deep', True):
            m = m.doit(**hints)
        if m.is_zero:
            return S.Zero
        if m.equals(S.One):
            return S.One
        if m.equals(S.Half):
            return exp(-pi)
        return self.func(m)

    def _eval_rewrite_as_elliptic_k(self, *args):
        m = self.args[0]
        return exp(-pi*elliptic_k(1 - m)/elliptic_k(m))

    _eval_rewrite_as_exp = _eval_rewrite_as_elliptic_k

    def fdiff(self, argindex=1):
        if argindex == 1:
            m = self.args[0]
            return pi**2*elliptic_nome_q(m)/(4*(m-1)*m*elliptic_k(m)**2)
        raise ArgumentIndexError(self, argindex)


def elliptic_nome_q_from_k(k):
    return elliptic_nome_q(k**2)

def elliptic_modulus_from_nome_q(q):
    return theta2(q)**2/theta3(q)**2
