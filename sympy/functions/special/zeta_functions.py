""" Riemann zeta and related function. """
from sympy.core import Function, S, C, sympify, pi
from sympy.core.function import ArgumentIndexError

###############################################################################
###################### LERCH TRANSCENDENT #####################################
###############################################################################

###############################################################################
###################### POLYLOGARITHM ##########################################
###############################################################################

###############################################################################
###################### HURWITZ GENERALIZED ZETA FUNCTION ######################
###############################################################################

class zeta(Function):
    r"""
    Hurwitz zeta function.

    For :math:`Re(a) > 0` and :math:`Re(s) > 1`, this function is defined as

    .. math:: \zeta(s, a) = \sum_{n=0}^\infty \frac{1}{(n + a)^s},

    where the standard choice of argument for :math:`n + a` is used. For fixed
    :math:`a` with :math:`Re(a) > 0` the Hurwitz zeta function admits a
    meromorphic continuation to all of :math:`\mathbb{C}`, it is an unbranched
    function with a simple pole at :math:`s = 1`.

    Analytic continuation to other :math:`a` is possible under some circumstances,
    but this is not typically done.

    **Examples**

    For :math:`a = 1` the Hurwitz zeta function reduces to the famous Riemann
    zeta function:

    .. math:: \zeta(s, 1) = \zeta(s) = \sum_{n=1}^\infty \frac{1}{n^s}.

    >>> from sympy import zeta
    >>> from sympy.abc import s
    >>> zeta(s, 1)
    zeta(s)

    The Riemann zeta function can also be expressed using the Dirichlet eta
    function:

    >>> from sympy import dirichlet_eta
    >>> zeta(s).rewrite(dirichlet_eta)
    dirichlet_eta(s)/(-2**(-s + 1) + 1)

    The Riemann zeta function at positive even integer and negative odd integer
    values is related to the Bernoulli numbers:

    >>> zeta(2)
    pi**2/6
    >>> zeta(4)
    pi**4/90
    >>> zeta(-1)
    -1/12

    At negative even integers the Riemann zeta function is zero:

    >>> zeta(-4)
    0

    No closed-form expressions are known at positive odd integers, but
    numerical evaluation is possible:

    >>> zeta(3).n()
    1.20205690315959

    The derivative of :math:`\zeta(s, a)` with respect to :math:`a` is easily
    computed:

    >>> from sympy.abc import a
    >>> zeta(s, a).diff(a)
    -s*zeta(s + 1, a)

    However the derivative with respect to :math:`s` has no useful closed form
    expression:

    >>> zeta(s, a).diff(s)
    Derivative(zeta(s, a), s)

    **See also:** :class:`dirichlet_eta`, :class:`lerchphi`, :class:`polylog`.

    **References**

    - http://dlmf.nist.gov/25.11
    - http://en.wikipedia.org/wiki/Hurwitz_zeta_function
    """

    nargs = (1, 2)

    @classmethod
    def eval(cls, z, a_=None):
        if a_ is None:
            z, a = map(sympify, (z, 1))
        else:
            z, a = map(sympify, (z, a_))

        if a.is_Number:
            if a is S.NaN:
                return S.NaN
            elif a is S.One and a_ is not None:
                return cls(z)
            # TODO Should a == 0 return S.NaN as well?

        if z.is_Number:
            if z is S.NaN:
                return S.NaN
            elif z is S.Infinity:
                return S.One
            elif z is S.Zero:
                if a.is_negative:
                    return S.Half - a - 1
                else:
                    return S.Half - a
            elif z is S.One:
                return S.ComplexInfinity
            elif z.is_Integer:
                if a.is_Integer:
                    if z.is_negative:
                        zeta = (-1)**z * C.bernoulli(-z+1)/(-z+1)
                    elif z.is_even:
                        B, F = C.bernoulli(z), C.factorial(z)
                        zeta = 2**(z-1) * abs(B) * pi**z / F
                    else:
                        return

                    if a.is_negative:
                        return zeta + C.harmonic(abs(a), z)
                    else:
                        return zeta - C.harmonic(a-1, z)

    def _eval_rewrite_as_dirichlet_eta(self, s, a=1):
        if a != 1:
            return self
        s = self.args[0]
        return dirichlet_eta(s)/(1 - 2**(1 - s))

    def fdiff(self, argindex=1):
        if len(self.args) == 2:
            s, a = self.args
        else:
            s, a = self.args + (1,)
        if argindex == 2:
            return -s*zeta(s + 1, a)
        else:
            raise ArgumentIndexError

class dirichlet_eta(Function):
    r"""
    Dirichlet eta function.

    For :math:`Re(s) > 0`, this function is defined as

    .. math:: \eta(s) = \sum_{n=1}^\infty \frac{(-1)^n}{n^s}.

    It admits a unique analytic continuation to all of :math:`\mathbb{C}`.
    It is an entire, unbranched function.

    **Examples**

    The Dirichlet eta function is closely related to the Riemann zeta function:

    >>> from sympy import dirichlet_eta, zeta
    >>> from sympy.abc import s
    >>> dirichlet_eta(s).rewrite(zeta)
    (-2**(-s + 1) + 1)*zeta(s)

    **See also:** :class:`zeta`

    **References**

    - http://en.wikipedia.org/wiki/Dirichlet_eta_function
    """
    nargs = 1

    @classmethod
    def eval(cls, s):
        if s == 1:
            return C.log(2)
        z = zeta(s)
        if not z.has(zeta):
            return (1-2**(1-s))*z

    def _eval_rewrite_as_zeta(self, s):
        return (1-2**(1-s)) * zeta(s)
