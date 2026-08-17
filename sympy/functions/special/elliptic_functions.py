"""Elliptic functions."""
from __future__ import annotations

from sympy.core import Add, S
from sympy.core.function import ArgumentIndexError, DefinedFunction
from sympy.core.numbers import pi
from sympy.functions.elementary.trigonometric import cos, sin


class jtheta(DefinedFunction):
    r"""
    The Jacobi theta function ``jtheta(n, z, q[, d])`` of index $n$.

    Explanation
    ===========

    The Jacobi theta functions are defined for $n = 1, 2, 3, 4$ by

    .. math::
        \vartheta_1(z, q) = 2 \sum_{k=0}^{\infty}
            (-1)^k q^{(k + 1/2)^2} \sin((2k + 1)z),

    .. math::
        \vartheta_2(z, q) = 2 \sum_{k=0}^{\infty}
            q^{(k + 1/2)^2} \cos((2k + 1)z),

    .. math::
        \vartheta_3(z, q) = 1 + 2 \sum_{k=1}^{\infty}
            q^{k^2} \cos(2kz),

    .. math::
        \vartheta_4(z, q) = 1 + 2 \sum_{k=1}^{\infty}
            (-1)^k q^{k^2} \cos(2kz).

    Here $z$ is the argument and $q$ is the nome, with $|q| < 1$. The
    optional fourth argument $d$ is a nonnegative integer denoting the order
    of differentiation with respect to $z$:

    .. math::
        \operatorname{jtheta}(n,z,q,d)
        = \frac{\partial^d}{\partial z^d}\vartheta_n(z,q).

    Omitting $d$ is equivalent to taking $d=0$.

    SymPy also supports differentiation with respect to the nome. Since
    $q=\exp(i\pi\tau)$, the theta-function heat equation [3]_ gives

    .. math::
        \frac{\partial}{\partial q}\vartheta_n(z,q)
        =-\frac{1}{4q}\frac{\partial^2}{\partial z^2}\vartheta_n(z,q).

    More generally, differentiating ``jtheta(n, z, q, d)`` with respect to
    $q$ returns ``-jtheta(n, z, q, d + 2)/(4*q)``.

    For integers $m$ and $n$, the functions are quasi-periodic on the
    period lattice [1]_:

    .. math::
        \begin{aligned}
        \vartheta_1(z+(m+n\tau)\pi\mid\tau)
            &=(-1)^{m+n}q^{-n^2}e^{-2inz}\vartheta_1(z\mid\tau),\\
        \vartheta_2(z+(m+n\tau)\pi\mid\tau)
            &=(-1)^m q^{-n^2}e^{-2inz}\vartheta_2(z\mid\tau),\\
        \vartheta_3(z+(m+n\tau)\pi\mid\tau)
            &=q^{-n^2}e^{-2inz}\vartheta_3(z\mid\tau),\\
        \vartheta_4(z+(m+n\tau)\pi\mid\tau)
            &=(-1)^n q^{-n^2}e^{-2inz}\vartheta_4(z\mid\tau).
        \end{aligned}

    Translation by half-periods interchanges the four functions. With
    $M=\exp(iz+i\pi\tau/4)$, the relations are

    .. math::
        \begin{aligned}
        \vartheta_1(z\mid\tau)
            &=-\vartheta_2(z+\tfrac{\pi}{2}\mid\tau)
             =-iM\vartheta_4(z+\tfrac{\pi\tau}{2}\mid\tau)
             =-iM\vartheta_3(z+\tfrac{\pi+\pi\tau}{2}\mid\tau),\\
        \vartheta_2(z\mid\tau)
            &=\vartheta_1(z+\tfrac{\pi}{2}\mid\tau)
             =M\vartheta_3(z+\tfrac{\pi\tau}{2}\mid\tau)
             =M\vartheta_4(z+\tfrac{\pi+\pi\tau}{2}\mid\tau),\\
        \vartheta_3(z\mid\tau)
            &=\vartheta_4(z+\tfrac{\pi}{2}\mid\tau)
             =M\vartheta_2(z+\tfrac{\pi\tau}{2}\mid\tau)
             =M\vartheta_1(z+\tfrac{\pi+\pi\tau}{2}\mid\tau),\\
        \vartheta_4(z\mid\tau)
            &=\vartheta_3(z+\tfrac{\pi}{2}\mid\tau)
             =-iM\vartheta_1(z+\tfrac{\pi\tau}{2}\mid\tau)
             =iM\vartheta_2(z+\tfrac{\pi+\pi\tau}{2}\mid\tau).
        \end{aligned}

    Examples
    ========

    >>> from sympy import diff, jtheta
    >>> from sympy.abc import n, q, z

    The first theta function is odd in its argument, whereas the other
    three are even:

    >>> jtheta(1, -z, q)
    -jtheta(1, z, q)
    >>> jtheta(3, -z, q)
    jtheta(3, z, q)

    Derivatives with respect to the argument are represented by an optional
    fourth argument:

    >>> diff(jtheta(n, z, q), z)
    jtheta(n, z, q, 1)
    >>> diff(jtheta(n, z, q), q)
    -jtheta(n, z, q, 2)/(4*q)

    The quasi-periodicity and half-period translation formulas can be
    checked numerically. SymPy does not currently apply these translations
    as symbolic simplifications:

    >>> from sympy import I, Rational, exp, pi
    >>> tau = 2*I/5
    >>> q = exp(I*pi*tau)
    >>> z = Rational(1, 3) + I/7
    >>> shifted = jtheta(4, z + pi*tau, q)
    >>> expected = -exp(-2*I*z)*jtheta(4, z, q)/q
    >>> abs((shifted - expected).evalf(30)) < 1e-25
    True
    >>> M = exp(I*z + I*pi*tau/4)
    >>> shifted = -I*M*jtheta(4, z + pi*tau/2, q)
    >>> abs((jtheta(1, z, q) - shifted).evalf(30)) < 1e-25
    True

    See Also
    ========

    elliptic_k

    References
    ==========

    .. [1] https://dlmf.nist.gov/20.2
    .. [2] https://mpmath.org/doc/current/functions/elliptic.html#jtheta
    .. [3] https://dlmf.nist.gov/20.13.E1

    """

    @classmethod
    def eval(cls, n, z, q, d=None):
        if n.is_Integer:
            if n not in (S.One, S(2), S(3), S(4)):
                raise ValueError("jtheta index must be 1, 2, 3 or 4")
        elif n.is_Float or n.is_integer is False:
            raise ValueError("jtheta index must be an integer")

        if d is not None:
            if d.is_Float or d.is_integer is False:
                raise ValueError("jtheta derivative order must be "
                                 "a nonnegative integer")
            if d.is_nonnegative is False:
                raise ValueError("jtheta derivative order must be "
                                 "a nonnegative integer")
            if d.is_zero:
                return cls(n, z, q)

        derivative_order = S.Zero if d is None else d

        if q.is_zero:
            if n in (S.One, S(2)):
                return S.Zero
            if n in (S(3), S(4)):
                if d is None:
                    return S.One
                if derivative_order.is_positive:
                    return S.Zero

        if z.is_zero:
            if n is S.One and derivative_order.is_even:
                return S.Zero
            if n in (S(2), S(3), S(4)) and derivative_order.is_odd:
                return S.Zero

        if z.could_extract_minus_sign():
            sign = None
            if n is S.One:
                if derivative_order.is_even:
                    sign = S.NegativeOne
                elif derivative_order.is_odd:
                    sign = S.One
            elif n in (S(2), S(3), S(4)):
                if derivative_order.is_even:
                    sign = S.One
                elif derivative_order.is_odd:
                    sign = S.NegativeOne
            if sign is not None:
                args = (n, -z, q) if d is None else (n, -z, q, d)
                return sign*cls(*args)

    @property
    def derivative_order(self):
        """Order of differentiation with respect to the argument."""
        if len(self.args) == 3:
            return S.Zero
        return self.args[3]

    def fdiff(self, argindex=2):
        n, z, q = self.args[:3]
        d = self.derivative_order
        if argindex == 2:
            return self.func(n, z, q, d + 1)
        if argindex == 3:
            return -self.func(n, z, q, d + 2)/(4*q)
        raise ArgumentIndexError(self, argindex)

    def _eval_nseries(self, x, n, logx, cdir=0):
        from sympy.series.order import Order

        index, z, q = self.args[:3]
        d = self.derivative_order

        if not q.has(x) or q.subs(x, S.Zero).is_zero is not True:
            return super()._eval_nseries(x, n, logx, cdir)

        if (index not in (S.One, S(2), S(3), S(4))
                or not d.is_Integer or d.is_nonnegative is not True):
            raise NotImplementedError(
                "jtheta series requires an explicit index and derivative "
                "order")

        q_leading = q.as_leading_term(x, logx=logx, cdir=cdir)
        _, q_order = q_leading.as_coeff_exponent(x)
        if q_order.is_number is not True or q_order.is_positive is not True:
            raise NotImplementedError(
                "jtheta series requires a nome with explicit positive order")

        terms = []
        if index in (S(3), S(4)):
            if d.is_zero:
                terms.append(S.One)
            k = 1
            while k**2*q_order < n:
                coefficient = 2*(2*k)**d*cos(2*k*z + d*pi/2)
                if index is S(4):
                    coefficient *= S.NegativeOne**k
                terms.append(coefficient*q**(k**2))
                k += 1
        else:
            k = 0
            exponent = (S(k) + S.Half)**2
            while exponent*q_order < n:
                harmonic = 2*k + 1
                if index is S.One:
                    coefficient = (2*S.NegativeOne**k*harmonic**d
                                   * sin(harmonic*z + d*pi/2))
                else:
                    coefficient = (2*harmonic**d
                                   * cos(harmonic*z + d*pi/2))
                terms.append(coefficient*q**exponent)
                k += 1
                exponent = (S(k) + S.Half)**2

        series = Add(*terms)
        return series._eval_nseries(x, n, logx, cdir) + Order(x**n, x)
