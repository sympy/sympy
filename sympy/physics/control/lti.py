from typing import Type

from sympy.core.add import Add
from sympy.core.basic import Basic
from sympy.core.containers import Tuple
from sympy.core.evalf import EvalfMixin
from sympy.core.expr import Expr
from sympy.core.function import expand
from sympy.core.logic import fuzzy_and
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.symbol import Dummy, Symbol
from sympy.core.sympify import sympify, _sympify
from sympy.matrices import ImmutableMatrix, eye
from sympy.matrices.expressions import MatMul, MatAdd
from sympy.polys import Poly, rootof
from sympy.polys.polyroots import roots
from sympy.polys.polytools import (cancel, degree)
from sympy.series import limit

from mpmath.libmp.libmpf import prec_to_dps

__all__ = ['TransferFunction', 'Series', 'MIMOSeries', 'Parallel', 'MIMOParallel',
    'Feedback', 'MIMOFeedback', 'TransferFunctionMatrix']


def _roots(poly, var):
    """ like roots, but works on higher-order polynomials. """
    r = roots(poly, var, multiple=True)
    n = degree(poly)
    if len(r) != n:
        r = [rootof(poly, var, k) for k in range(n)]
    return r


class LinearTimeInvariant(Basic, EvalfMixin):
    """A common class for all the Linear Time-Invariant Dynamical Systems."""

    _clstype: Type

    # Users should not directly interact with this class.
    def __new__(cls, *system, **kwargs):
        if cls is LinearTimeInvariant:
            raise NotImplementedError('The LTICommon class is not meant to be used directly.')
        return super(LinearTimeInvariant, cls).__new__(cls, *system, **kwargs)

    @classmethod
    def _check_args(cls, args):
        if not args:
            raise ValueError("Atleast 1 argument must be passed.")
        if not all(isinstance(arg, cls._clstype) for arg in args):
            raise TypeError(f"All arguments must be of type {cls._clstype}.")
        var_set = {arg.var for arg in args}
        if len(var_set) != 1:
            raise ValueError("All transfer functions should use the same complex variable"
                f" of the Laplace transform. {len(var_set)} different values found.")

    @property
    def is_SISO(self):
        """Returns `True` if the passed LTI system is SISO else returns False."""
        return self._is_SISO


class SISOLinearTimeInvariant(LinearTimeInvariant):
    """A common class for all the SISO Linear Time-Invariant Dynamical Systems."""
    # Users should not directly interact with this class.
    _is_SISO = True


class MIMOLinearTimeInvariant(LinearTimeInvariant):
    """A common class for all the MIMO Linear Time-Invariant Dynamical Systems."""
    # Users should not directly interact with this class.
    _is_SISO = False


SISOLinearTimeInvariant._clstype = SISOLinearTimeInvariant
MIMOLinearTimeInvariant._clstype = MIMOLinearTimeInvariant


def _check_other_SISO(func):
    def wrapper(*args, **kwargs):
        if not isinstance(args[-1], SISOLinearTimeInvariant):
            return NotImplemented
        else:
            return func(*args, **kwargs)
    return wrapper


def _check_other_MIMO(func):
    def wrapper(*args, **kwargs):
        if not isinstance(args[-1], MIMOLinearTimeInvariant):
            return NotImplemented
        else:
            return func(*args, **kwargs)
    return wrapper


class TransferFunction(SISOLinearTimeInvariant):
    r"""
    A class for representing LTI (Linear, time-invariant) systems that can be strictly described
    by ratio of polynomials in the Laplace transform complex variable. The arguments
    are ``num``, ``den``, and ``var``, where ``num`` and ``den`` are numerator and
    denominator polynomials of the ``TransferFunction`` respectively, and the third argument is
    a complex variable of the Laplace transform used by these polynomials of the transfer function.
    ``num`` and ``den`` can be either polynomials or numbers, whereas ``var``
    has to be a :py:class:`~.Symbol`.

    Explanation
    ===========

    Generally, a dynamical system representing a physical model can be described in terms of Linear
    Ordinary Differential Equations like -

            $\small{b_{m}y^{\left(m\right)}+b_{m-1}y^{\left(m-1\right)}+\dots+b_{1}y^{\left(1\right)}+b_{0}y=
            a_{n}x^{\left(n\right)}+a_{n-1}x^{\left(n-1\right)}+\dots+a_{1}x^{\left(1\right)}+a_{0}x}$

    Here, $x$ is the input signal and $y$ is the output signal and superscript on both is the order of derivative
    (not exponent). Derivative is taken with respect to the independent variable, $t$. Also, generally $m$ is greater
    than $n$.

    It is not feasible to analyse the properties of such systems in their native form therefore, we use
    mathematical tools like Laplace transform to get a better perspective. Taking the Laplace transform
    of both the sides in the equation (at zero initial conditions), we get -

            $\small{\mathcal{L}[b_{m}y^{\left(m\right)}+b_{m-1}y^{\left(m-1\right)}+\dots+b_{1}y^{\left(1\right)}+b_{0}y]=
            \mathcal{L}[a_{n}x^{\left(n\right)}+a_{n-1}x^{\left(n-1\right)}+\dots+a_{1}x^{\left(1\right)}+a_{0}x]}$

    Using the linearity property of Laplace transform and also considering zero initial conditions
    (i.e. $\small{y(0^{-}) = 0}$, $\small{y'(0^{-}) = 0}$ and so on), the equation
    above gets translated to -

            $\small{b_{m}\mathcal{L}[y^{\left(m\right)}]+\dots+b_{1}\mathcal{L}[y^{\left(1\right)}]+b_{0}\mathcal{L}[y]=
            a_{n}\mathcal{L}[x^{\left(n\right)}]+\dots+a_{1}\mathcal{L}[x^{\left(1\right)}]+a_{0}\mathcal{L}[x]}$

    Now, applying Derivative property of Laplace transform,

            $\small{b_{m}s^{m}\mathcal{L}[y]+\dots+b_{1}s\mathcal{L}[y]+b_{0}\mathcal{L}[y]=
            a_{n}s^{n}\mathcal{L}[x]+\dots+a_{1}s\mathcal{L}[x]+a_{0}\mathcal{L}[x]}$

    Here, the superscript on $s$ is **exponent**. Note that the zero initial conditions assumption, mentioned above, is very important
    and cannot be ignored otherwise the dynamical system cannot be considered time-independent and the simplified equation above
    cannot be reached.

    Collecting $\mathcal{L}[y]$ and $\mathcal{L}[x]$ terms from both the sides and taking the ratio
    $\frac{ \mathcal{L}\left\{y\right\} }{ \mathcal{L}\left\{x\right\} }$, we get the typical rational form of transfer
    function.

    The numerator of the transfer function is, therefore, the Laplace transform of the output signal
    (The signals are represented as functions of time) and similarly, the denominator
    of the transfer function is the Laplace transform of the input signal. It is also a convention
    to denote the input and output signal's Laplace transform with capital alphabets like shown below.

            $H(s) = \frac{Y(s)}{X(s)} = \frac{ \mathcal{L}\left\{y(t)\right\} }{ \mathcal{L}\left\{x(t)\right\} }$

    $s$, also known as complex frequency, is a complex variable in the Laplace domain. It corresponds to the
    equivalent variable $t$, in the time domain. Transfer functions are sometimes also referred to as the Laplace
    transform of the system's impulse response. Transfer function, $H$, is represented as a rational
    function in $s$ like,

            $H(s) =\ \frac{a_{n}s^{n}+a_{n-1}s^{n-1}+\dots+a_{1}s+a_{0}}{b_{m}s^{m}+b_{m-1}s^{m-1}+\dots+b_{1}s+b_{0}}$

    Parameters
    ==========

    num : Expr, Number
        The numerator polynomial of the transfer function.
    den : Expr, Number
        The denominator polynomial of the transfer function.
    var : Symbol
        Complex variable of the Laplace transform used by the
        polynomials of the transfer function.

    Raises
    ======

    TypeError
        When ``var`` is not a Symbol or when ``num`` or ``den`` is not a
        number or a polynomial.
    ValueError
        When ``den`` is zero.

    Examples
    ========

    >>> from sympy.abc import s, p, a
    >>> from sympy.physics.control.lti import TransferFunction
    >>> tf1 = TransferFunction(s + a, s**2 + s + 1, s)
    >>> tf1
    TransferFunction(a + s, s**2 + s + 1, s)
    >>> tf1.num
    a + s
    >>> tf1.den
    s**2 + s + 1
    >>> tf1.var
    s
    >>> tf1.args
    (a + s, s**2 + s + 1, s)

    Any complex variable can be used for ``var``.

    >>> tf2 = TransferFunction(a*p**3 - a*p**2 + s*p, p + a**2, p)
    >>> tf2
    TransferFunction(a*p**3 - a*p**2 + p*s, a**2 + p, p)
    >>> tf3 = TransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p)
    >>> tf3
    TransferFunction((p - 1)*(p + 3), (p - 1)*(p + 5), p)

    To negate a transfer function the ``-`` operator can be prepended:

    >>> tf4 = TransferFunction(-a + s, p**2 + s, p)
    >>> -tf4
    TransferFunction(a - s, p**2 + s, p)
    >>> tf5 = TransferFunction(s**4 - 2*s**3 + 5*s + 4, s + 4, s)
    >>> -tf5
    TransferFunction(-s**4 + 2*s**3 - 5*s - 4, s + 4, s)

    You can use a float or an integer (or other constants) as numerator and denominator:

    >>> tf6 = TransferFunction(1/2, 4, s)
    >>> tf6.num
    0.500000000000000
    >>> tf6.den
    4
    >>> tf6.var
    s
    >>> tf6.args
    (0.5, 4, s)

    You can take the integer power of a transfer function using the ``**`` operator:

    >>> tf7 = TransferFunction(s + a, s - a, s)
    >>> tf7**3
    TransferFunction((a + s)**3, (-a + s)**3, s)
    >>> tf7**0
    TransferFunction(1, 1, s)
    >>> tf8 = TransferFunction(p + 4, p - 3, p)
    >>> tf8**-1
    TransferFunction(p - 3, p + 4, p)

    Addition, subtraction, and multiplication of transfer functions can form
    unevaluated ``Series`` or ``Parallel`` objects.

    >>> tf9 = TransferFunction(s + 1, s**2 + s + 1, s)
    >>> tf10 = TransferFunction(s - p, s + 3, s)
    >>> tf11 = TransferFunction(4*s**2 + 2*s - 4, s - 1, s)
    >>> tf12 = TransferFunction(1 - s, s**2 + 4, s)
    >>> tf9 + tf10
    Parallel(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(-p + s, s + 3, s))
    >>> tf10 - tf11
    Parallel(TransferFunction(-p + s, s + 3, s), TransferFunction(-4*s**2 - 2*s + 4, s - 1, s))
    >>> tf9 * tf10
    Series(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(-p + s, s + 3, s))
    >>> tf10 - (tf9 + tf12)
    Parallel(TransferFunction(-p + s, s + 3, s), TransferFunction(-s - 1, s**2 + s + 1, s), TransferFunction(s - 1, s**2 + 4, s))
    >>> tf10 - (tf9 * tf12)
    Parallel(TransferFunction(-p + s, s + 3, s), Series(TransferFunction(-1, 1, s), TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(1 - s, s**2 + 4, s)))
    >>> tf11 * tf10 * tf9
    Series(TransferFunction(4*s**2 + 2*s - 4, s - 1, s), TransferFunction(-p + s, s + 3, s), TransferFunction(s + 1, s**2 + s + 1, s))
    >>> tf9 * tf11 + tf10 * tf12
    Parallel(Series(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(4*s**2 + 2*s - 4, s - 1, s)), Series(TransferFunction(-p + s, s + 3, s), TransferFunction(1 - s, s**2 + 4, s)))
    >>> (tf9 + tf12) * (tf10 + tf11)
    Series(Parallel(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(1 - s, s**2 + 4, s)), Parallel(TransferFunction(-p + s, s + 3, s), TransferFunction(4*s**2 + 2*s - 4, s - 1, s)))

    These unevaluated ``Series`` or ``Parallel`` objects can convert into the
    resultant transfer function using ``.doit()`` method or by ``.rewrite(TransferFunction)``.

    >>> ((tf9 + tf10) * tf12).doit()
    TransferFunction((1 - s)*((-p + s)*(s**2 + s + 1) + (s + 1)*(s + 3)), (s + 3)*(s**2 + 4)*(s**2 + s + 1), s)
    >>> (tf9 * tf10 - tf11 * tf12).rewrite(TransferFunction)
    TransferFunction(-(1 - s)*(s + 3)*(s**2 + s + 1)*(4*s**2 + 2*s - 4) + (-p + s)*(s - 1)*(s + 1)*(s**2 + 4), (s - 1)*(s + 3)*(s**2 + 4)*(s**2 + s + 1), s)

    See Also
    ========

    Feedback, Series, Parallel

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Transfer_function
    .. [2] https://en.wikipedia.org/wiki/Laplace_transform

    """
    def __new__(cls, num, den, var):
        num, den = _sympify(num), _sympify(den)

        if not isinstance(var, Symbol):
            raise TypeError("Variable input must be a Symbol.")

        if den == 0:
            raise ValueError("TransferFunction cannot have a zero denominator.")

        if (((isinstance(num, Expr) and num.has(Symbol)) or num.is_number) and
            ((isinstance(den, Expr) and den.has(Symbol)) or den.is_number)):
            obj = super(TransferFunction, cls).__new__(cls, num, den, var)
            obj._num = num
            obj._den = den
            obj._var = var
            return obj

        else:
            raise TypeError("Unsupported type for numerator or denominator of TransferFunction.")

    @classmethod
    def from_rational_expression(cls, expr, var=None):
        r"""
        Creates a new ``TransferFunction`` efficiently from a rational expression.

        Parameters
        ==========

        expr : Expr, Number
            The rational expression representing the ``TransferFunction``.
        var : Symbol, optional
            Complex variable of the Laplace transform used by the
            polynomials of the transfer function.

        Raises
        ======

        ValueError
            When ``expr`` is of type ``Number`` and optional parameter ``var``
            is not passed.

            When ``expr`` has more than one variables and an optional parameter
            ``var`` is not passed.
        ZeroDivisionError
            When denominator of ``expr`` is zero or it has ``ComplexInfinity``
            in its numerator.

        Examples
        ========

        >>> from sympy.abc import s, p, a
        >>> from sympy.physics.control.lti import TransferFunction
        >>> expr1 = (s + 5)/(3*s**2 + 2*s + 1)
        >>> tf1 = TransferFunction.from_rational_expression(expr1)
        >>> tf1
        TransferFunction(s + 5, 3*s**2 + 2*s + 1, s)
        >>> expr2 = (a*p**3 - a*p**2 + s*p)/(p + a**2)  # Expr with more than one variables
        >>> tf2 = TransferFunction.from_rational_expression(expr2, p)
        >>> tf2
        TransferFunction(a*p**3 - a*p**2 + p*s, a**2 + p, p)

        In case of conflict between two or more variables in a expression, SymPy will
        raise a ``ValueError``, if ``var`` is not passed by the user.

        >>> tf = TransferFunction.from_rational_expression((a + a*s)/(s**2 + s + 1))
        Traceback (most recent call last):
        ...
        ValueError: Conflicting values found for positional argument `var` ({a, s}). Specify it manually.

        This can be corrected by specifying the ``var`` parameter manually.

        >>> tf = TransferFunction.from_rational_expression((a + a*s)/(s**2 + s + 1), s)
        >>> tf
        TransferFunction(a*s + a, s**2 + s + 1, s)

        ``var`` also need to be specified when ``expr`` is a ``Number``

        >>> tf3 = TransferFunction.from_rational_expression(10, s)
        >>> tf3
        TransferFunction(10, 1, s)

        """
        expr = _sympify(expr)
        if var is None:
            _free_symbols = expr.free_symbols
            _len_free_symbols = len(_free_symbols)
            if _len_free_symbols == 1:
                var = list(_free_symbols)[0]
            elif _len_free_symbols == 0:
                raise ValueError("Positional argument `var` not found in the TransferFunction defined. Specify it manually.")
            else:
                raise ValueError("Conflicting values found for positional argument `var` ({}). Specify it manually.".format(_free_symbols))

        _num, _den = expr.as_numer_denom()
        if _den == 0 or _num.has(S.ComplexInfinity):
            raise ZeroDivisionError("TransferFunction cannot have a zero denominator.")
        return cls(_num, _den, var)

    @property
    def num(self):
        """
        Returns the numerator polynomial of the transfer function.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction
        >>> G1 = TransferFunction(s**2 + p*s + 3, s - 4, s)
        >>> G1.num
        p*s + s**2 + 3
        >>> G2 = TransferFunction((p + 5)*(p - 3), (p - 3)*(p + 1), p)
        >>> G2.num
        (p - 3)*(p + 5)

        """
        return self._num

    @property
    def den(self):
        """
        Returns the denominator polynomial of the transfer function.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction
        >>> G1 = TransferFunction(s + 4, p**3 - 2*p + 4, s)
        >>> G1.den
        p**3 - 2*p + 4
        >>> G2 = TransferFunction(3, 4, s)
        >>> G2.den
        4

        """
        return self._den

    @property
    def var(self):
        """
        Returns the complex variable of the Laplace transform used by the polynomials of
        the transfer function.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction
        >>> G1 = TransferFunction(p**2 + 2*p + 4, p - 6, p)
        >>> G1.var
        p
        >>> G2 = TransferFunction(0, s - 5, s)
        >>> G2.var
        s

        """
        return self._var

    def _eval_subs(self, old, new):
        arg_num = self.num.subs(old, new)
        arg_den = self.den.subs(old, new)
        argnew = TransferFunction(arg_num, arg_den, self.var)
        return self if old == self.var else argnew

    def _eval_evalf(self, prec):
        return TransferFunction(
            self.num._eval_evalf(prec),
            self.den._eval_evalf(prec),
            self.var)

    def _eval_simplify(self, **kwargs):
        tf = cancel(Mul(self.num, 1/self.den, evaluate=False), expand=False).as_numer_denom()
        num_, den_ = tf[0], tf[1]
        return TransferFunction(num_, den_, self.var)

    def expand(self):
        """
        Returns the transfer function with numerator and denominator
        in expanded form.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction
        >>> G1 = TransferFunction((a - s)**2, (s**2 + a)**2, s)
        >>> G1.expand()
        TransferFunction(a**2 - 2*a*s + s**2, a**2 + 2*a*s**2 + s**4, s)
        >>> G2 = TransferFunction((p + 3*b)*(p - b), (p - b)*(p + 2*b), p)
        >>> G2.expand()
        TransferFunction(-3*b**2 + 2*b*p + p**2, -2*b**2 + b*p + p**2, p)

        """
        return TransferFunction(expand(self.num), expand(self.den), self.var)

    def dc_gain(self):
        """
        Computes the gain of the response as the frequency approaches zero.

        The DC gain is infinite for systems with pure integrators.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction(s + 3, s**2 - 9, s)
        >>> tf1.dc_gain()
        -1/3
        >>> tf2 = TransferFunction(p**2, p - 3 + p**3, p)
        >>> tf2.dc_gain()
        0
        >>> tf3 = TransferFunction(a*p**2 - b, s + b, s)
        >>> tf3.dc_gain()
        (a*p**2 - b)/b
        >>> tf4 = TransferFunction(1, s, s)
        >>> tf4.dc_gain()
        oo

        """
        m = Mul(self.num, Pow(self.den, -1, evaluate=False), evaluate=False)
        return limit(m, self.var, 0)

    def poles(self):
        """
        Returns the poles of a transfer function.

        Examples
        ========

        >>> from sympy.abc import s, p, a
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p)
        >>> tf1.poles()
        [-5, 1]
        >>> tf2 = TransferFunction((1 - s)**2, (s**2 + 1)**2, s)
        >>> tf2.poles()
        [I, I, -I, -I]
        >>> tf3 = TransferFunction(s**2, a*s + p, s)
        >>> tf3.poles()
        [-p/a]

        """
        return _roots(Poly(self.den, self.var), self.var)

    def zeros(self):
        """
        Returns the zeros of a transfer function.

        Examples
        ========

        >>> from sympy.abc import s, p, a
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p)
        >>> tf1.zeros()
        [-3, 1]
        >>> tf2 = TransferFunction((1 - s)**2, (s**2 + 1)**2, s)
        >>> tf2.zeros()
        [1, 1]
        >>> tf3 = TransferFunction(s**2, a*s + p, s)
        >>> tf3.zeros()
        [0, 0]

        """
        return _roots(Poly(self.num, self.var), self.var)

    def is_stable(self):
        """
        Returns True if the transfer function is asymptotically stable; else False.

        This would not check the marginal or conditional stability of the system.

        Examples
        ========

        >>> from sympy.abc import s, p, a
        >>> from sympy import symbols
        >>> from sympy.physics.control.lti import TransferFunction
        >>> q, r = symbols('q, r', negative=True)
        >>> tf1 = TransferFunction((1 - s)**2, (s + 1)**2, s)
        >>> tf1.is_stable()
        True
        >>> tf2 = TransferFunction((1 - p)**2, (s**2 + 1)**2, s)
        >>> tf2.is_stable()
        False
        >>> tf3 = TransferFunction(4, q*s - r, s)
        >>> tf3.is_stable()
        False
        >>> tf4 = TransferFunction(p + 1, a*p - s**2, p)
        >>> tf4.is_stable() is None   # Not enough info about the symbols to determine stability
        True

        """
        return fuzzy_and(pole.as_real_imag()[0].is_negative for pole in self.poles())

    def __add__(self, other):
        if isinstance(other, (TransferFunction, Series)):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            return Parallel(self, other)
        elif isinstance(other, Parallel):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            arg_list = list(other.args)
            return Parallel(self, *arg_list)
        else:
            raise ValueError("TransferFunction cannot be added with {}.".
                format(type(other)))

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if isinstance(other, (TransferFunction, Series)):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            return Parallel(self, -other)
        elif isinstance(other, Parallel):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            arg_list = [-i for i in list(other.args)]
            return Parallel(self, *arg_list)
        else:
            raise ValueError("{} cannot be subtracted from a TransferFunction."
                .format(type(other)))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if isinstance(other, (TransferFunction, Parallel)):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            return Series(self, other)
        elif isinstance(other, Series):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            arg_list = list(other.args)
            return Series(self, *arg_list)
        else:
            raise ValueError("TransferFunction cannot be multiplied with {}."
                .format(type(other)))

    __rmul__ = __mul__

    def __truediv__(self, other):
        if (isinstance(other, Parallel) and len(other.args) == 2 and isinstance(other.args[0], TransferFunction)
            and isinstance(other.args[1], (Series, TransferFunction))):

            if not self.var == other.var:
                raise ValueError("Both TransferFunction and Parallel should use the"
                    " same complex variable of the Laplace transform.")
            if other.args[1] == self:
                # plant and controller with unit feedback.
                return Feedback(self, other.args[0])
            other_arg_list = list(other.args[1].args) if isinstance(other.args[1], Series) else other.args[1]
            if other_arg_list == other.args[1]:
                return Feedback(self, other_arg_list)
            elif self in other_arg_list:
                other_arg_list.remove(self)
            else:
                return Feedback(self, Series(*other_arg_list))

            if len(other_arg_list) == 1:
                return Feedback(self, *other_arg_list)
            else:
                return Feedback(self, Series(*other_arg_list))
        else:
            raise ValueError("TransferFunction cannot be divided by {}.".
                format(type(other)))

    __rtruediv__ = __truediv__

    def __pow__(self, p):
        p = sympify(p)
        if not p.is_Integer:
            raise ValueError("Exponent must be an integer.")
        if p is S.Zero:
            return TransferFunction(1, 1, self.var)
        elif p > 0:
            num_, den_ = self.num**p, self.den**p
        else:
            p = abs(p)
            num_, den_ = self.den**p, self.num**p

        return TransferFunction(num_, den_, self.var)

    def __neg__(self):
        return TransferFunction(-self.num, self.den, self.var)

    @property
    def is_proper(self):
        """
        Returns True if degree of the numerator polynomial is less than
        or equal to degree of the denominator polynomial, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
        >>> tf1.is_proper
        False
        >>> tf2 = TransferFunction(p**2 - 4*p, p**3 + 3*p + 2, p)
        >>> tf2.is_proper
        True

        """
        return degree(self.num, self.var) <= degree(self.den, self.var)

    @property
    def is_strictly_proper(self):
        """
        Returns True if degree of the numerator polynomial is strictly less
        than degree of the denominator polynomial, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf1.is_strictly_proper
        False
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> tf2.is_strictly_proper
        True

        """
        return degree(self.num, self.var) < degree(self.den, self.var)

    @property
    def is_biproper(self):
        """
        Returns True if degree of the numerator polynomial is equal to
        degree of the denominator polynomial, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf1.is_biproper
        True
        >>> tf2 = TransferFunction(p**2, p + a, p)
        >>> tf2.is_biproper
        False

        """
        return degree(self.num, self.var) == degree(self.den, self.var)

    def to_expr(self):
        """
        Converts a ``TransferFunction`` object to SymPy Expr.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction
        >>> from sympy import Expr
        >>> tf1 = TransferFunction(s, a*s**2 + 1, s)
        >>> tf1.to_expr()
        s/(a*s**2 + 1)
        >>> isinstance(_, Expr)
        True
        >>> tf2 = TransferFunction(1, (p + 3*b)*(b - p), p)
        >>> tf2.to_expr()
        1/((b - p)*(3*b + p))
        >>> tf3 = TransferFunction((s - 2)*(s - 3), (s - 1)*(s - 2)*(s - 3), s)
        >>> tf3.to_expr()
        ((s - 3)*(s - 2))/(((s - 3)*(s - 2)*(s - 1)))

        """

        if self.num != 1:
            return Mul(self.num, Pow(self.den, -1, evaluate=False), evaluate=False)
        else:
            return Pow(self.den, -1, evaluate=False)


def _flatten_args(args, _cls):
    temp_args = []
    for arg in args:
        if isinstance(arg, _cls):
            temp_args.extend(arg.args)
        else:
            temp_args.append(arg)
    return tuple(temp_args)


def _dummify_args(_arg, var):
    dummy_dict = {}
    dummy_arg_list = []

    for arg in _arg:
        _s = Dummy()
        dummy_dict[_s] = var
        dummy_arg = arg.subs({var: _s})
        dummy_arg_list.append(dummy_arg)

    return dummy_arg_list, dummy_dict


class Series(SISOLinearTimeInvariant):
    r"""
    A class for representing a series configuration of SISO systems.

    Parameters
    ==========

    args : SISOLinearTimeInvariant
        SISO systems in a series configuration.
    evaluate : Boolean, Keyword
        When passed ``True``, returns the equivalent
        ``Series(*args).doit()``. Set to ``False`` by default.

    Raises
    ======

    ValueError
        When no argument is passed.

        ``var`` attribute is not same for every system.
    TypeError
        Any of the passed ``*args`` has unsupported type

        A combination of SISO and MIMO systems is
        passed. There should be homogeneity in the
        type of systems passed, SISO in this case.

    Examples
    ========

    >>> from sympy.abc import s, p, a, b
    >>> from sympy.physics.control.lti import TransferFunction, Series, Parallel
    >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
    >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
    >>> tf3 = TransferFunction(p**2, p + s, s)
    >>> S1 = Series(tf1, tf2)
    >>> S1
    Series(TransferFunction(a*p**2 + b*s, -p + s, s), TransferFunction(s**3 - 2, s**4 + 5*s + 6, s))
    >>> S1.var
    s
    >>> S2 = Series(tf2, Parallel(tf3, -tf1))
    >>> S2
    Series(TransferFunction(s**3 - 2, s**4 + 5*s + 6, s), Parallel(TransferFunction(p**2, p + s, s), TransferFunction(-a*p**2 - b*s, -p + s, s)))
    >>> S2.var
    s
    >>> S3 = Series(Parallel(tf1, tf2), Parallel(tf2, tf3))
    >>> S3
    Series(Parallel(TransferFunction(a*p**2 + b*s, -p + s, s), TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)), Parallel(TransferFunction(s**3 - 2, s**4 + 5*s + 6, s), TransferFunction(p**2, p + s, s)))
    >>> S3.var
    s

    You can get the resultant transfer function by using ``.doit()`` method:

    >>> S3 = Series(tf1, tf2, -tf3)
    >>> S3.doit()
    TransferFunction(-p**2*(s**3 - 2)*(a*p**2 + b*s), (-p + s)*(p + s)*(s**4 + 5*s + 6), s)
    >>> S4 = Series(tf2, Parallel(tf1, -tf3))
    >>> S4.doit()
    TransferFunction((s**3 - 2)*(-p**2*(-p + s) + (p + s)*(a*p**2 + b*s)), (-p + s)*(p + s)*(s**4 + 5*s + 6), s)

    Notes
    =====

    All the transfer functions should use the same complex variable
    ``var`` of the Laplace transform.

    See Also
    ========

    MIMOSeries, Parallel, TransferFunction, Feedback

    """
    def __new__(cls, *args, evaluate=False):

        args = _flatten_args(args, Series)
        cls._check_args(args)
        obj = super().__new__(cls, *args)

        return obj.doit() if evaluate else obj

    @property
    def var(self):
        """
        Returns the complex variable used by all the transfer functions.

        Examples
        ========

        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, Series, Parallel
        >>> G1 = TransferFunction(p**2 + 2*p + 4, p - 6, p)
        >>> G2 = TransferFunction(p, 4 - p, p)
        >>> G3 = TransferFunction(0, p**4 - 1, p)
        >>> Series(G1, G2).var
        p
        >>> Series(-G3, Parallel(G1, G2)).var
        p

        """
        return self.args[0].var

    def doit(self, **hints):
        """
        Returns the resultant transfer function obtained after evaluating
        the transfer functions in series configuration.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Series
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> Series(tf2, tf1).doit()
        TransferFunction((s**3 - 2)*(a*p**2 + b*s), (-p + s)*(s**4 + 5*s + 6), s)
        >>> Series(-tf1, -tf2).doit()
        TransferFunction((2 - s**3)*(-a*p**2 - b*s), (-p + s)*(s**4 + 5*s + 6), s)

        """

        _num_arg = (arg.doit().num for arg in self.args)
        _den_arg = (arg.doit().den for arg in self.args)
        res_num = Mul(*_num_arg, evaluate=True)
        res_den = Mul(*_den_arg, evaluate=True)
        return TransferFunction(res_num, res_den, self.var)

    def _eval_rewrite_as_TransferFunction(self, *args, **kwargs):
        return self.doit()

    @_check_other_SISO
    def __add__(self, other):

        if isinstance(other, Parallel):
            arg_list = list(other.args)
            return Parallel(self, *arg_list)

        return Parallel(self, other)

    __radd__ = __add__

    @_check_other_SISO
    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    @_check_other_SISO
    def __mul__(self, other):

        arg_list = list(self.args)
        return Series(*arg_list, other)

    def __truediv__(self, other):
        if (isinstance(other, Parallel) and len(other.args) == 2
            and isinstance(other.args[0], TransferFunction) and isinstance(other.args[1], Series)):

            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            self_arg_list = set(list(self.args))
            other_arg_list = set(list(other.args[1].args))
            res = list(self_arg_list ^ other_arg_list)
            if len(res) == 0:
                return Feedback(self, other.args[0])
            elif len(res) == 1:
                return Feedback(self, *res)
            else:
                return Feedback(self, Series(*res))
        else:
            raise ValueError("This transfer function expression is invalid.")

    def __neg__(self):
        return Series(TransferFunction(-1, 1, self.var), self)

    def to_expr(self):
        """Returns the equivalent ``Expr`` object."""
        return Mul(*(arg.to_expr() for arg in self.args), evaluate=False)

    @property
    def is_proper(self):
        """
        Returns True if degree of the numerator polynomial of the resultant transfer
        function is less than or equal to degree of the denominator polynomial of
        the same, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Series
        >>> tf1 = TransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
        >>> tf2 = TransferFunction(p**2 - 4*p, p**3 + 3*s + 2, s)
        >>> tf3 = TransferFunction(s, s**2 + s + 1, s)
        >>> S1 = Series(-tf2, tf1)
        >>> S1.is_proper
        False
        >>> S2 = Series(tf1, tf2, tf3)
        >>> S2.is_proper
        True

        """
        return self.doit().is_proper

    @property
    def is_strictly_proper(self):
        """
        Returns True if degree of the numerator polynomial of the resultant transfer
        function is strictly less than degree of the denominator polynomial of
        the same, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Series
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**2 + 5*s + 6, s)
        >>> tf3 = TransferFunction(1, s**2 + s + 1, s)
        >>> S1 = Series(tf1, tf2)
        >>> S1.is_strictly_proper
        False
        >>> S2 = Series(tf1, tf2, tf3)
        >>> S2.is_strictly_proper
        True

        """
        return self.doit().is_strictly_proper

    @property
    def is_biproper(self):
        r"""
        Returns True if degree of the numerator polynomial of the resultant transfer
        function is equal to degree of the denominator polynomial of
        the same, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Series
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(p, s**2, s)
        >>> tf3 = TransferFunction(s**2, 1, s)
        >>> S1 = Series(tf1, -tf2)
        >>> S1.is_biproper
        False
        >>> S2 = Series(tf2, tf3)
        >>> S2.is_biproper
        True

        """
        return self.doit().is_biproper


def _mat_mul_compatible(*args):
    """To check whether shapes are compatible for matrix mul."""
    return all(args[i].num_outputs == args[i+1].num_inputs for i in range(len(args)-1))


class MIMOSeries(MIMOLinearTimeInvariant):
    r"""
    A class for representing a series configuration of MIMO systems.

    Parameters
    ==========

    args : MIMOLinearTimeInvariant
        MIMO systems in a series configuration.
    evaluate : Boolean, Keyword
        When passed ``True``, returns the equivalent
        ``MIMOSeries(*args).doit()``. Set to ``False`` by default.

    Raises
    ======

    ValueError
        When no argument is passed.

        ``var`` attribute is not same for every system.

        ``num_outputs`` of the MIMO system is not equal to the
        ``num_inputs`` of its adjacent MIMO system. (Matrix
        multiplication constraint, basically)
    TypeError
        Any of the passed ``*args`` has unsupported type

        A combination of SISO and MIMO systems is
        passed. There should be homogeneity in the
        type of systems passed, MIMO in this case.

    Examples
    ========

    >>> from sympy.abc import s
    >>> from sympy.physics.control.lti import MIMOSeries, TransferFunctionMatrix
    >>> from sympy import Matrix, pprint
    >>> mat_a = Matrix([[5*s], [5]])  # 2 Outputs 1 Input
    >>> mat_b = Matrix([[5, 1/(6*s**2)]])  # 1 Output 2 Inputs
    >>> mat_c = Matrix([[1, s], [5/s, 1]])  # 2 Outputs 2 Inputs
    >>> tfm_a = TransferFunctionMatrix.from_Matrix(mat_a, s)
    >>> tfm_b = TransferFunctionMatrix.from_Matrix(mat_b, s)
    >>> tfm_c = TransferFunctionMatrix.from_Matrix(mat_c, s)
    >>> MIMOSeries(tfm_c, tfm_b, tfm_a)
    MIMOSeries(TransferFunctionMatrix(((TransferFunction(1, 1, s), TransferFunction(s, 1, s)), (TransferFunction(5, s, s), TransferFunction(1, 1, s)))), TransferFunctionMatrix(((TransferFunction(5, 1, s), TransferFunction(1, 6*s**2, s)),)), TransferFunctionMatrix(((TransferFunction(5*s, 1, s),), (TransferFunction(5, 1, s),))))
    >>> pprint(_, use_unicode=False)  #  For Better Visualization
    [5*s]                 [1  s]
    [---]    [5   1  ]    [-  -]
    [ 1 ]    [-  ----]    [1  1]
    [   ]   *[1     2]   *[    ]
    [ 5 ]    [   6*s ]{t} [5  1]
    [ - ]                 [-  -]
    [ 1 ]{t}              [s  1]{t}
    >>> MIMOSeries(tfm_c, tfm_b, tfm_a).doit()
    TransferFunctionMatrix(((TransferFunction(150*s**4 + 25*s, 6*s**3, s), TransferFunction(150*s**4 + 5*s, 6*s**2, s)), (TransferFunction(150*s**3 + 25, 6*s**3, s), TransferFunction(150*s**3 + 5, 6*s**2, s))))
    >>> pprint(_, use_unicode=False)  # (2 Inputs -A-> 2 Outputs) -> (2 Inputs -B-> 1 Output) -> (1 Input -C-> 2 Outputs) is equivalent to (2 Inputs -Series Equivalent-> 2 Outputs).
    [     4              4      ]
    [150*s  + 25*s  150*s  + 5*s]
    [-------------  ------------]
    [        3             2    ]
    [     6*s           6*s     ]
    [                           ]
    [      3              3     ]
    [ 150*s  + 25    150*s  + 5 ]
    [ -----------    ---------- ]
    [        3             2    ]
    [     6*s           6*s     ]{t}

    Notes
    =====

    All the transfer function matrices should use the same complex variable ``var`` of the Laplace transform.

    ``MIMOSeries(A, B)`` is not equivalent to ``A*B``. It is always in the reverse order, that is ``B*A``.

    See Also
    ========

    Series, MIMOParallel

    """
    def __new__(cls, *args, evaluate=False):

        cls._check_args(args)

        if _mat_mul_compatible(*args):
            obj = super().__new__(cls, *args)

        else:
            raise ValueError("Number of input signals do not match the number"
                " of output signals of adjacent systems for some args.")

        return obj.doit() if evaluate else obj

    @property
    def var(self):
        """
        Returns the complex variable used by all the transfer functions.

        Examples
        ========

        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, MIMOSeries, TransferFunctionMatrix
        >>> G1 = TransferFunction(p**2 + 2*p + 4, p - 6, p)
        >>> G2 = TransferFunction(p, 4 - p, p)
        >>> G3 = TransferFunction(0, p**4 - 1, p)
        >>> tfm_1 = TransferFunctionMatrix([[G1, G2, G3]])
        >>> tfm_2 = TransferFunctionMatrix([[G1], [G2], [G3]])
        >>> MIMOSeries(tfm_2, tfm_1).var
        p

        """
        return self.args[0].var

    @property
    def num_inputs(self):
        """Returns the number of input signals of the series system."""
        return self.args[0].num_inputs

    @property
    def num_outputs(self):
        """Returns the number of output signals of the series system."""
        return self.args[-1].num_outputs

    @property
    def shape(self):
        """Returns the shape of the equivalent MIMO system."""
        return self.num_outputs, self.num_inputs

    def doit(self, cancel=False, **kwargs):
        """
        Returns the resultant transfer function matrix obtained after evaluating
        the MIMO systems arranged in a series configuration.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, MIMOSeries, TransferFunctionMatrix
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> tfm1 = TransferFunctionMatrix([[tf1, tf2], [tf2, tf2]])
        >>> tfm2 = TransferFunctionMatrix([[tf2, tf1], [tf1, tf1]])
        >>> MIMOSeries(tfm2, tfm1).doit()
        TransferFunctionMatrix(((TransferFunction(2*(-p + s)*(s**3 - 2)*(a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)**2*(s**4 + 5*s + 6)**2, s), TransferFunction((-p + s)**2*(s**3 - 2)*(a*p**2 + b*s) + (-p + s)*(a*p**2 + b*s)**2*(s**4 + 5*s + 6), (-p + s)**3*(s**4 + 5*s + 6), s)), (TransferFunction((-p + s)*(s**3 - 2)**2*(s**4 + 5*s + 6) + (s**3 - 2)*(a*p**2 + b*s)*(s**4 + 5*s + 6)**2, (-p + s)*(s**4 + 5*s + 6)**3, s), TransferFunction(2*(s**3 - 2)*(a*p**2 + b*s), (-p + s)*(s**4 + 5*s + 6), s))))

        """
        _arg = (arg.doit()._expr_mat for arg in reversed(self.args))

        if cancel:
            res = MatMul(*_arg, evaluate=True)
            return TransferFunctionMatrix.from_Matrix(res, self.var)

        _dummy_args, _dummy_dict = _dummify_args(_arg, self.var)
        res = MatMul(*_dummy_args, evaluate=True)
        temp_tfm = TransferFunctionMatrix.from_Matrix(res, self.var)
        return temp_tfm.subs(_dummy_dict)

    def _eval_rewrite_as_TransferFunctionMatrix(self, *args, **kwargs):
        return self.doit()

    @_check_other_MIMO
    def __add__(self, other):

        if isinstance(other, MIMOParallel):
            arg_list = list(other.args)
            return MIMOParallel(self, *arg_list)

        return MIMOParallel(self, other)

    __radd__ = __add__

    @_check_other_MIMO
    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    @_check_other_MIMO
    def __mul__(self, other):

        if isinstance(other, MIMOSeries):
            self_arg_list = list(self.args)
            other_arg_list = list(other.args)
            return MIMOSeries(*other_arg_list, *self_arg_list)  # A*B = MIMOSeries(B, A)

        arg_list = list(self.args)
        return MIMOSeries(other, *arg_list)

    def __neg__(self):
        arg_list = list(self.args)
        arg_list[0] = -arg_list[0]
        return MIMOSeries(*arg_list)


class Parallel(SISOLinearTimeInvariant):
    r"""
    A class for representing a parallel configuration of SISO systems.

    Parameters
    ==========

    args : SISOLinearTimeInvariant
        SISO systems in a parallel arrangement.
    evaluate : Boolean, Keyword
        When passed ``True``, returns the equivalent
        ``Parallel(*args).doit()``. Set to ``False`` by default.

    Raises
    ======

    ValueError
        When no argument is passed.

        ``var`` attribute is not same for every system.
    TypeError
        Any of the passed ``*args`` has unsupported type

        A combination of SISO and MIMO systems is
        passed. There should be homogeneity in the
        type of systems passed.

    Examples
    ========

    >>> from sympy.abc import s, p, a, b
    >>> from sympy.physics.control.lti import TransferFunction, Parallel, Series
    >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
    >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
    >>> tf3 = TransferFunction(p**2, p + s, s)
    >>> P1 = Parallel(tf1, tf2)
    >>> P1
    Parallel(TransferFunction(a*p**2 + b*s, -p + s, s), TransferFunction(s**3 - 2, s**4 + 5*s + 6, s))
    >>> P1.var
    s
    >>> P2 = Parallel(tf2, Series(tf3, -tf1))
    >>> P2
    Parallel(TransferFunction(s**3 - 2, s**4 + 5*s + 6, s), Series(TransferFunction(p**2, p + s, s), TransferFunction(-a*p**2 - b*s, -p + s, s)))
    >>> P2.var
    s
    >>> P3 = Parallel(Series(tf1, tf2), Series(tf2, tf3))
    >>> P3
    Parallel(Series(TransferFunction(a*p**2 + b*s, -p + s, s), TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)), Series(TransferFunction(s**3 - 2, s**4 + 5*s + 6, s), TransferFunction(p**2, p + s, s)))
    >>> P3.var
    s

    You can get the resultant transfer function by using ``.doit()`` method:

    >>> Parallel(tf1, tf2, -tf3).doit()
    TransferFunction(-p**2*(-p + s)*(s**4 + 5*s + 6) + (-p + s)*(p + s)*(s**3 - 2) + (p + s)*(a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)*(p + s)*(s**4 + 5*s + 6), s)
    >>> Parallel(tf2, Series(tf1, -tf3)).doit()
    TransferFunction(-p**2*(a*p**2 + b*s)*(s**4 + 5*s + 6) + (-p + s)*(p + s)*(s**3 - 2), (-p + s)*(p + s)*(s**4 + 5*s + 6), s)

    Notes
    =====

    All the transfer functions should use the same complex variable
    ``var`` of the Laplace transform.

    See Also
    ========

    Series, TransferFunction, Feedback

    """
    def __new__(cls, *args, evaluate=False):

        args = _flatten_args(args, Parallel)
        cls._check_args(args)
        obj = super().__new__(cls, *args)

        return obj.doit() if evaluate else obj

    @property
    def var(self):
        """
        Returns the complex variable used by all the transfer functions.

        Examples
        ========

        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, Parallel, Series
        >>> G1 = TransferFunction(p**2 + 2*p + 4, p - 6, p)
        >>> G2 = TransferFunction(p, 4 - p, p)
        >>> G3 = TransferFunction(0, p**4 - 1, p)
        >>> Parallel(G1, G2).var
        p
        >>> Parallel(-G3, Series(G1, G2)).var
        p

        """
        return self.args[0].var

    def doit(self, **hints):
        """
        Returns the resultant transfer function obtained after evaluating
        the transfer functions in parallel configuration.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Parallel
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> Parallel(tf2, tf1).doit()
        TransferFunction((-p + s)*(s**3 - 2) + (a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s)
        >>> Parallel(-tf1, -tf2).doit()
        TransferFunction((2 - s**3)*(-p + s) + (-a*p**2 - b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s)

        """

        _arg = (arg.doit().to_expr() for arg in self.args)
        res = Add(*_arg).as_numer_denom()
        return TransferFunction(*res, self.var)

    def _eval_rewrite_as_TransferFunction(self, *args, **kwargs):
        return self.doit()

    @_check_other_SISO
    def __add__(self, other):

        self_arg_list = list(self.args)
        return Parallel(*self_arg_list, other)

    __radd__ = __add__

    @_check_other_SISO
    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    @_check_other_SISO
    def __mul__(self, other):

        if isinstance(other, Series):
            arg_list = list(other.args)
            return Series(self, *arg_list)

        return Series(self, other)

    def __neg__(self):
        return Series(TransferFunction(-1, 1, self.var), self)

    def to_expr(self):
        """Returns the equivalent ``Expr`` object."""
        return Add(*(arg.to_expr() for arg in self.args), evaluate=False)

    @property
    def is_proper(self):
        """
        Returns True if degree of the numerator polynomial of the resultant transfer
        function is less than or equal to degree of the denominator polynomial of
        the same, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Parallel
        >>> tf1 = TransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
        >>> tf2 = TransferFunction(p**2 - 4*p, p**3 + 3*s + 2, s)
        >>> tf3 = TransferFunction(s, s**2 + s + 1, s)
        >>> P1 = Parallel(-tf2, tf1)
        >>> P1.is_proper
        False
        >>> P2 = Parallel(tf2, tf3)
        >>> P2.is_proper
        True

        """
        return self.doit().is_proper

    @property
    def is_strictly_proper(self):
        """
        Returns True if degree of the numerator polynomial of the resultant transfer
        function is strictly less than degree of the denominator polynomial of
        the same, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Parallel
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> tf3 = TransferFunction(s, s**2 + s + 1, s)
        >>> P1 = Parallel(tf1, tf2)
        >>> P1.is_strictly_proper
        False
        >>> P2 = Parallel(tf2, tf3)
        >>> P2.is_strictly_proper
        True

        """
        return self.doit().is_strictly_proper

    @property
    def is_biproper(self):
        """
        Returns True if degree of the numerator polynomial of the resultant transfer
        function is equal to degree of the denominator polynomial of
        the same, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Parallel
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(p**2, p + s, s)
        >>> tf3 = TransferFunction(s, s**2 + s + 1, s)
        >>> P1 = Parallel(tf1, -tf2)
        >>> P1.is_biproper
        True
        >>> P2 = Parallel(tf2, tf3)
        >>> P2.is_biproper
        False

        """
        return self.doit().is_biproper


class MIMOParallel(MIMOLinearTimeInvariant):
    r"""
    A class for representing a parallel configuration of MIMO systems.

    Parameters
    ==========

    args : MIMOLinearTimeInvariant
        MIMO Systems in a parallel arrangement.
    evaluate : Boolean, Keyword
        When passed ``True``, returns the equivalent
        ``MIMOParallel(*args).doit()``. Set to ``False`` by default.

    Raises
    ======

    ValueError
        When no argument is passed.

        ``var`` attribute is not same for every system.

        All MIMO systems passed do not have same shape.
    TypeError
        Any of the passed ``*args`` has unsupported type

        A combination of SISO and MIMO systems is
        passed. There should be homogeneity in the
        type of systems passed, MIMO in this case.

    Examples
    ========

    >>> from sympy.abc import s
    >>> from sympy.physics.control.lti import TransferFunctionMatrix, MIMOParallel
    >>> from sympy import Matrix, pprint
    >>> expr_1 = 1/s
    >>> expr_2 = s/(s**2-1)
    >>> expr_3 = (2 + s)/(s**2 - 1)
    >>> expr_4 = 5
    >>> tfm_a = TransferFunctionMatrix.from_Matrix(Matrix([[expr_1, expr_2], [expr_3, expr_4]]), s)
    >>> tfm_b = TransferFunctionMatrix.from_Matrix(Matrix([[expr_2, expr_1], [expr_4, expr_3]]), s)
    >>> tfm_c = TransferFunctionMatrix.from_Matrix(Matrix([[expr_3, expr_4], [expr_1, expr_2]]), s)
    >>> MIMOParallel(tfm_a, tfm_b, tfm_c)
    MIMOParallel(TransferFunctionMatrix(((TransferFunction(1, s, s), TransferFunction(s, s**2 - 1, s)), (TransferFunction(s + 2, s**2 - 1, s), TransferFunction(5, 1, s)))), TransferFunctionMatrix(((TransferFunction(s, s**2 - 1, s), TransferFunction(1, s, s)), (TransferFunction(5, 1, s), TransferFunction(s + 2, s**2 - 1, s)))), TransferFunctionMatrix(((TransferFunction(s + 2, s**2 - 1, s), TransferFunction(5, 1, s)), (TransferFunction(1, s, s), TransferFunction(s, s**2 - 1, s)))))
    >>> pprint(_, use_unicode=False)  #  For Better Visualization
    [  1       s   ]      [  s       1   ]      [s + 2     5   ]
    [  -     ------]      [------    -   ]      [------    -   ]
    [  s      2    ]      [ 2        s   ]      [ 2        1   ]
    [        s  - 1]      [s  - 1        ]      [s  - 1        ]
    [              ]    + [              ]    + [              ]
    [s + 2     5   ]      [  5     s + 2 ]      [  1       s   ]
    [------    -   ]      [  -     ------]      [  -     ------]
    [ 2        1   ]      [  1      2    ]      [  s      2    ]
    [s  - 1        ]{t}   [        s  - 1]{t}   [        s  - 1]{t}
    >>> MIMOParallel(tfm_a, tfm_b, tfm_c).doit()
    TransferFunctionMatrix(((TransferFunction(s**2 + s*(2*s + 2) - 1, s*(s**2 - 1), s), TransferFunction(2*s**2 + 5*s*(s**2 - 1) - 1, s*(s**2 - 1), s)), (TransferFunction(s**2 + s*(s + 2) + 5*s*(s**2 - 1) - 1, s*(s**2 - 1), s), TransferFunction(5*s**2 + 2*s - 3, s**2 - 1, s))))
    >>> pprint(_, use_unicode=False)
    [       2                              2       / 2    \    ]
    [      s  + s*(2*s + 2) - 1         2*s  + 5*s*\s  - 1/ - 1]
    [      --------------------         -----------------------]
    [             / 2    \                       / 2    \      ]
    [           s*\s  - 1/                     s*\s  - 1/      ]
    [                                                          ]
    [ 2                   / 2    \             2               ]
    [s  + s*(s + 2) + 5*s*\s  - 1/ - 1      5*s  + 2*s - 3     ]
    [---------------------------------      --------------     ]
    [              / 2    \                      2             ]
    [            s*\s  - 1/                     s  - 1         ]{t}

    Notes
    =====

    All the transfer function matrices should use the same complex variable
    ``var`` of the Laplace transform.

    See Also
    ========

    Parallel, MIMOSeries

    """
    def __new__(cls, *args, evaluate=False):

        args = _flatten_args(args, MIMOParallel)

        cls._check_args(args)

        if any(arg.shape != args[0].shape for arg in args):
             raise TypeError("Shape of all the args is not equal.")

        obj = super().__new__(cls, *args)

        return obj.doit() if evaluate else obj

    @property
    def var(self):
        """
        Returns the complex variable used by all the systems.

        Examples
        ========

        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, MIMOParallel
        >>> G1 = TransferFunction(p**2 + 2*p + 4, p - 6, p)
        >>> G2 = TransferFunction(p, 4 - p, p)
        >>> G3 = TransferFunction(0, p**4 - 1, p)
        >>> G4 = TransferFunction(p**2, p**2 - 1, p)
        >>> tfm_a = TransferFunctionMatrix([[G1, G2], [G3, G4]])
        >>> tfm_b = TransferFunctionMatrix([[G2, G1], [G4, G3]])
        >>> MIMOParallel(tfm_a, tfm_b).var
        p

        """
        return self.args[0].var

    @property
    def num_inputs(self):
        """Returns the number of input signals of the parallel system."""
        return self.args[0].num_inputs

    @property
    def num_outputs(self):
        """Returns the number of output signals of the parallel system."""
        return self.args[0].num_outputs

    @property
    def shape(self):
        """Returns the shape of the equivalent MIMO system."""
        return self.num_outputs, self.num_inputs

    def doit(self, **hints):
        """
        Returns the resultant transfer function matrix obtained after evaluating
        the MIMO systems arranged in a parallel configuration.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, MIMOParallel, TransferFunctionMatrix
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> tfm_1 = TransferFunctionMatrix([[tf1, tf2], [tf2, tf1]])
        >>> tfm_2 = TransferFunctionMatrix([[tf2, tf1], [tf1, tf2]])
        >>> MIMOParallel(tfm_1, tfm_2).doit()
        TransferFunctionMatrix(((TransferFunction((-p + s)*(s**3 - 2) + (a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s), TransferFunction((-p + s)*(s**3 - 2) + (a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s)), (TransferFunction((-p + s)*(s**3 - 2) + (a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s), TransferFunction((-p + s)*(s**3 - 2) + (a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s))))

        """
        _arg = (arg.doit()._expr_mat for arg in self.args)
        res = MatAdd(*_arg, evaluate=True)
        return TransferFunctionMatrix.from_Matrix(res, self.var)

    def _eval_rewrite_as_TransferFunctionMatrix(self, *args, **kwargs):
        return self.doit()

    @_check_other_MIMO
    def __add__(self, other):

        self_arg_list = list(self.args)
        return MIMOParallel(*self_arg_list, other)

    __radd__ = __add__

    @_check_other_MIMO
    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    @_check_other_MIMO
    def __mul__(self, other):

        if isinstance(other, MIMOSeries):
            arg_list = list(other.args)
            return MIMOSeries(*arg_list, self)

        return MIMOSeries(other, self)

    def __neg__(self):
        arg_list = [-arg for arg in list(self.args)]
        return MIMOParallel(*arg_list)


class Feedback(SISOLinearTimeInvariant):
    r"""
    A class for representing closed-loop feedback interconnection between two
    SISO input/output systems.

    The first argument, ``sys1``, is the feedforward part of the closed-loop
    system or in simple words, the dynamical model representing the process
    to be controlled. The second argument, ``sys2``, is the feedback system
    and controls the fed back signal to ``sys1``. Both ``sys1`` and ``sys2``
    can either be ``Series`` or ``TransferFunction`` objects.

    Parameters
    ==========

    sys1 : Series, TransferFunction
        The feedforward path system.
    sys2 : Series, TransferFunction, optional
        The feedback path system (often a feedback controller).
        It is the model sitting on the feedback path.

        If not specified explicitly, the sys2 is
        assumed to be unit (1.0) transfer function.
    sign : int, optional
        The sign of feedback. Can either be ``1``
        (for positive feedback) or ``-1`` (for negative feedback).
        Default value is `-1`.

    Raises
    ======

    ValueError
        When ``sys1`` and ``sys2`` are not using the
        same complex variable of the Laplace transform.

        When a combination of ``sys1`` and ``sys2`` yields
        zero denominator.

    TypeError
        When either ``sys1`` or ``sys2`` is not a ``Series`` or a
        ``TransferFunction`` object.

    Examples
    ========

    >>> from sympy.abc import s
    >>> from sympy.physics.control.lti import TransferFunction, Feedback
    >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
    >>> controller = TransferFunction(5*s - 10, s + 7, s)
    >>> F1 = Feedback(plant, controller)
    >>> F1
    Feedback(TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s), TransferFunction(5*s - 10, s + 7, s), -1)
    >>> F1.var
    s
    >>> F1.args
    (TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s), TransferFunction(5*s - 10, s + 7, s), -1)

    You can get the feedforward and feedback path systems by using ``.sys1`` and ``.sys2`` respectively.

    >>> F1.sys1
    TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
    >>> F1.sys2
    TransferFunction(5*s - 10, s + 7, s)

    You can get the resultant closed loop transfer function obtained by negative feedback
    interconnection using ``.doit()`` method.

    >>> F1.doit()
    TransferFunction((s + 7)*(s**2 - 4*s + 2)*(3*s**2 + 7*s - 3), ((s + 7)*(s**2 - 4*s + 2) + (5*s - 10)*(3*s**2 + 7*s - 3))*(s**2 - 4*s + 2), s)
    >>> G = TransferFunction(2*s**2 + 5*s + 1, s**2 + 2*s + 3, s)
    >>> C = TransferFunction(5*s + 10, s + 10, s)
    >>> F2 = Feedback(G*C, TransferFunction(1, 1, s))
    >>> F2.doit()
    TransferFunction((s + 10)*(5*s + 10)*(s**2 + 2*s + 3)*(2*s**2 + 5*s + 1), (s + 10)*((s + 10)*(s**2 + 2*s + 3) + (5*s + 10)*(2*s**2 + 5*s + 1))*(s**2 + 2*s + 3), s)

    To negate a ``Feedback`` object, the ``-`` operator can be prepended:

    >>> -F1
    Feedback(TransferFunction(-3*s**2 - 7*s + 3, s**2 - 4*s + 2, s), TransferFunction(10 - 5*s, s + 7, s), -1)
    >>> -F2
    Feedback(Series(TransferFunction(-1, 1, s), TransferFunction(2*s**2 + 5*s + 1, s**2 + 2*s + 3, s), TransferFunction(5*s + 10, s + 10, s)), TransferFunction(-1, 1, s), -1)

    See Also
    ========

    MIMOFeedback, Series, Parallel

    """
    def __new__(cls, sys1, sys2=None, sign=-1):
        if not sys2:
            sys2 = TransferFunction(1, 1, sys1.var)

        if not (isinstance(sys1, (TransferFunction, Series))
            and isinstance(sys2, (TransferFunction, Series))):
            raise TypeError("Unsupported type for `sys1` or `sys2` of Feedback.")

        if sign not in [-1, 1]:
            raise ValueError("Unsupported type for feedback. `sign` arg should "
                "either be 1 (positive feedback loop) or -1 (negative feedback loop).")

        if Mul(sys1.to_expr(), sys2.to_expr()).simplify() == sign:
            raise ValueError("The equivalent system will have zero denominator.")

        if sys1.var != sys2.var:
            raise ValueError("Both `sys1` and `sys2` should be using the"
                " same complex variable.")

        return super().__new__(cls, sys1, sys2, _sympify(sign))

    @property
    def sys1(self):
        """
        Returns the feedforward system of the feedback interconnection.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> controller = TransferFunction(5*s - 10, s + 7, s)
        >>> F1 = Feedback(plant, controller)
        >>> F1.sys1
        TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> G = TransferFunction(2*s**2 + 5*s + 1, p**2 + 2*p + 3, p)
        >>> C = TransferFunction(5*p + 10, p + 10, p)
        >>> P = TransferFunction(1 - s, p + 2, p)
        >>> F2 = Feedback(TransferFunction(1, 1, p), G*C*P)
        >>> F2.sys1
        TransferFunction(1, 1, p)

        """
        return self.args[0]

    @property
    def sys2(self):
        """
        Returns the feedback controller of the feedback interconnection.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> controller = TransferFunction(5*s - 10, s + 7, s)
        >>> F1 = Feedback(plant, controller)
        >>> F1.sys2
        TransferFunction(5*s - 10, s + 7, s)
        >>> G = TransferFunction(2*s**2 + 5*s + 1, p**2 + 2*p + 3, p)
        >>> C = TransferFunction(5*p + 10, p + 10, p)
        >>> P = TransferFunction(1 - s, p + 2, p)
        >>> F2 = Feedback(TransferFunction(1, 1, p), G*C*P)
        >>> F2.sys2
        Series(TransferFunction(2*s**2 + 5*s + 1, p**2 + 2*p + 3, p), TransferFunction(5*p + 10, p + 10, p), TransferFunction(1 - s, p + 2, p))

        """
        return self.args[1]

    @property
    def var(self):
        """
        Returns the complex variable of the Laplace transform used by all
        the transfer functions involved in the feedback interconnection.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> controller = TransferFunction(5*s - 10, s + 7, s)
        >>> F1 = Feedback(plant, controller)
        >>> F1.var
        s
        >>> G = TransferFunction(2*s**2 + 5*s + 1, p**2 + 2*p + 3, p)
        >>> C = TransferFunction(5*p + 10, p + 10, p)
        >>> P = TransferFunction(1 - s, p + 2, p)
        >>> F2 = Feedback(TransferFunction(1, 1, p), G*C*P)
        >>> F2.var
        p

        """
        return self.sys1.var

    @property
    def sign(self):
        """
        Returns the type of MIMO Feedback model. ``1``
        for Positive and ``-1`` for Negative.
        """
        return self.args[2]

    @property
    def sensitivity(self):
        """
        Returns the sensitivity function of the feedback loop.

        Sensitivity of a Feedback system is the ratio
        of change in the open loop gain to the change in
        the closed loop gain.

        .. note::
            This method would not return the complementary
            sensitivity function.

        Examples
        ========

        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> C = TransferFunction(5*p + 10, p + 10, p)
        >>> P = TransferFunction(1 - p, p + 2, p)
        >>> F_1 = Feedback(P, C)
        >>> F_1.sensitivity
        1/((1 - p)*(5*p + 10)/((p + 2)*(p + 10)) + 1)

        """

        return 1/(1 - self.sign*self.sys1.to_expr()*self.sys2.to_expr())

    def doit(self, cancel=False, expand=False, **hints):
        """
        Returns the resultant transfer function obtained by the
        feedback interconnection.

        Examples
        ========

        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> controller = TransferFunction(5*s - 10, s + 7, s)
        >>> F1 = Feedback(plant, controller)
        >>> F1.doit()
        TransferFunction((s + 7)*(s**2 - 4*s + 2)*(3*s**2 + 7*s - 3), ((s + 7)*(s**2 - 4*s + 2) + (5*s - 10)*(3*s**2 + 7*s - 3))*(s**2 - 4*s + 2), s)
        >>> G = TransferFunction(2*s**2 + 5*s + 1, s**2 + 2*s + 3, s)
        >>> F2 = Feedback(G, TransferFunction(1, 1, s))
        >>> F2.doit()
        TransferFunction((s**2 + 2*s + 3)*(2*s**2 + 5*s + 1), (s**2 + 2*s + 3)*(3*s**2 + 7*s + 4), s)

        Use kwarg ``expand=True`` to expand the resultant transfer function.
        Use ``cancel=True`` to cancel out the common terms in numerator and
        denominator.

        >>> F2.doit(cancel=True, expand=True)
        TransferFunction(2*s**2 + 5*s + 1, 3*s**2 + 7*s + 4, s)
        >>> F2.doit(expand=True)
        TransferFunction(2*s**4 + 9*s**3 + 17*s**2 + 17*s + 3, 3*s**4 + 13*s**3 + 27*s**2 + 29*s + 12, s)

        """
        arg_list = list(self.sys1.args) if isinstance(self.sys1, Series) else [self.sys1]
        # F_n and F_d are resultant TFs of num and den of Feedback.
        F_n, unit = self.sys1.doit(), TransferFunction(1, 1, self.sys1.var)
        if self.sign == -1:
            F_d = Parallel(unit, Series(self.sys2, *arg_list)).doit()
        else:
            F_d = Parallel(unit, -Series(self.sys2, *arg_list)).doit()

        _resultant_tf = TransferFunction(F_n.num * F_d.den, F_n.den * F_d.num, F_n.var)

        if cancel:
            _resultant_tf = _resultant_tf.simplify()

        if expand:
            _resultant_tf = _resultant_tf.expand()

        return _resultant_tf

    def _eval_rewrite_as_TransferFunction(self, num, den, sign, **kwargs):
        return self.doit()

    def __neg__(self):
        return Feedback(-self.sys1, -self.sys2, self.sign)


def _is_invertible(a, b, sign):
    """
    Checks whether a given pair of MIMO
    systems passed is invertible or not.
    """
    _mat = eye(a.num_outputs) - sign*(a.doit()._expr_mat)*(b.doit()._expr_mat)
    _det = _mat.det()

    return _det != 0


class MIMOFeedback(MIMOLinearTimeInvariant):
    r"""
    A class for representing closed-loop feedback interconnection between two
    MIMO input/output systems.

    Parameters
    ==========

    sys1 : MIMOSeries, TransferFunctionMatrix
        The MIMO system placed on the feedforward path.
    sys2 : MIMOSeries, TransferFunctionMatrix
        The system placed on the feedback path
        (often a feedback controller).
    sign : int, optional
        The sign of feedback. Can either be ``1``
        (for positive feedback) or ``-1`` (for negative feedback).
        Default value is `-1`.

    Raises
    ======

    ValueError
        When ``sys1`` and ``sys2`` are not using the
        same complex variable of the Laplace transform.

        Forward path model should have an equal number of inputs/outputs
        to the feedback path outputs/inputs.

        When product of ``sys1`` and ``sys2`` is not a square matrix.

        When the equivalent MIMO system is not invertible.

    TypeError
        When either ``sys1`` or ``sys2`` is not a ``MIMOSeries`` or a
        ``TransferFunctionMatrix`` object.

    Examples
    ========

    >>> from sympy import Matrix, pprint
    >>> from sympy.abc import s
    >>> from sympy.physics.control.lti import TransferFunctionMatrix, MIMOFeedback
    >>> plant_mat = Matrix([[1, 1/s], [0, 1]])
    >>> controller_mat = Matrix([[10, 0], [0, 10]])  # Constant Gain
    >>> plant = TransferFunctionMatrix.from_Matrix(plant_mat, s)
    >>> controller = TransferFunctionMatrix.from_Matrix(controller_mat, s)
    >>> feedback = MIMOFeedback(plant, controller)  # Negative Feedback (default)
    >>> pprint(feedback, use_unicode=False)
    /    [1  1]    [10  0 ]   \-1   [1  1]
    |    [-  -]    [--  - ]   |     [-  -]
    |    [1  s]    [1   1 ]   |     [1  s]
    |I + [    ]   *[      ]   |   * [    ]
    |    [0  1]    [0   10]   |     [0  1]
    |    [-  -]    [-   --]   |     [-  -]
    \    [1  1]{t} [1   1 ]{t}/     [1  1]{t}

    To get the equivalent system matrix, use either ``doit`` or ``rewrite`` method.

    >>> pprint(feedback.doit(), use_unicode=False)
    [1     1  ]
    [--  -----]
    [11  121*s]
    [         ]
    [0    1   ]
    [-    --  ]
    [1    11  ]{t}

    To negate the ``MIMOFeedback`` object, use ``-`` operator.

    >>> neg_feedback = -feedback
    >>> pprint(neg_feedback.doit(), use_unicode=False)
    [-1    -1  ]
    [---  -----]
    [ 11  121*s]
    [          ]
    [ 0    -1  ]
    [ -    --- ]
    [ 1     11 ]{t}

    See Also
    ========

    Feedback, MIMOSeries, MIMOParallel

    """
    def __new__(cls, sys1, sys2, sign=-1):
        if not (isinstance(sys1, (TransferFunctionMatrix, MIMOSeries))
            and isinstance(sys2, (TransferFunctionMatrix, MIMOSeries))):
            raise TypeError("Unsupported type for `sys1` or `sys2` of MIMO Feedback.")

        if sys1.num_inputs != sys2.num_outputs or \
            sys1.num_outputs != sys2.num_inputs:
            raise ValueError("Product of `sys1` and `sys2` "
                "must yield a square matrix.")

        if sign not in (-1, 1):
            raise ValueError("Unsupported type for feedback. `sign` arg should "
                "either be 1 (positive feedback loop) or -1 (negative feedback loop).")

        if not _is_invertible(sys1, sys2, sign):
            raise ValueError("Non-Invertible system inputted.")
        if sys1.var != sys2.var:
            raise ValueError("Both `sys1` and `sys2` should be using the"
                " same complex variable.")

        return super().__new__(cls, sys1, sys2, _sympify(sign))

    @property
    def sys1(self):
        r"""
        Returns the system placed on the feedforward path of the MIMO feedback interconnection.

        Examples
        ========

        >>> from sympy import pprint
        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, MIMOFeedback
        >>> tf1 = TransferFunction(s**2 + s + 1, s**2 - s + 1, s)
        >>> tf2 = TransferFunction(1, s, s)
        >>> tf3 = TransferFunction(1, 1, s)
        >>> sys1 = TransferFunctionMatrix([[tf1, tf2], [tf2, tf1]])
        >>> sys2 = TransferFunctionMatrix([[tf3, tf3], [tf3, tf2]])
        >>> F_1 = MIMOFeedback(sys1, sys2, 1)
        >>> F_1.sys1
        TransferFunctionMatrix(((TransferFunction(s**2 + s + 1, s**2 - s + 1, s), TransferFunction(1, s, s)), (TransferFunction(1, s, s), TransferFunction(s**2 + s + 1, s**2 - s + 1, s))))
        >>> pprint(_, use_unicode=False)
        [ 2                    ]
        [s  + s + 1      1     ]
        [----------      -     ]
        [ 2              s     ]
        [s  - s + 1            ]
        [                      ]
        [             2        ]
        [    1       s  + s + 1]
        [    -       ----------]
        [    s        2        ]
        [            s  - s + 1]{t}

        """
        return self.args[0]

    @property
    def sys2(self):
        r"""
        Returns the feedback controller of the MIMO feedback interconnection.

        Examples
        ========

        >>> from sympy import pprint
        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, MIMOFeedback
        >>> tf1 = TransferFunction(s**2, s**3 - s + 1, s)
        >>> tf2 = TransferFunction(1, s, s)
        >>> tf3 = TransferFunction(1, 1, s)
        >>> sys1 = TransferFunctionMatrix([[tf1, tf2], [tf2, tf1]])
        >>> sys2 = TransferFunctionMatrix([[tf1, tf3], [tf3, tf2]])
        >>> F_1 = MIMOFeedback(sys1, sys2)
        >>> F_1.sys2
        TransferFunctionMatrix(((TransferFunction(s**2, s**3 - s + 1, s), TransferFunction(1, 1, s)), (TransferFunction(1, 1, s), TransferFunction(1, s, s))))
        >>> pprint(_, use_unicode=False)
        [     2       ]
        [    s       1]
        [----------  -]
        [ 3          1]
        [s  - s + 1   ]
        [             ]
        [    1       1]
        [    -       -]
        [    1       s]{t}

        """
        return self.args[1]

    @property
    def var(self):
        r"""
        Returns the complex variable of the Laplace transform used by all
        the transfer functions involved in the MIMO feedback loop.

        Examples
        ========

        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, MIMOFeedback
        >>> tf1 = TransferFunction(p, 1 - p, p)
        >>> tf2 = TransferFunction(1, p, p)
        >>> tf3 = TransferFunction(1, 1, p)
        >>> sys1 = TransferFunctionMatrix([[tf1, tf2], [tf2, tf1]])
        >>> sys2 = TransferFunctionMatrix([[tf1, tf3], [tf3, tf2]])
        >>> F_1 = MIMOFeedback(sys1, sys2, 1)  # Positive feedback
        >>> F_1.var
        p

        """
        return self.sys1.var

    @property
    def sign(self):
        r"""
        Returns the type of feedback interconnection of two models. ``1``
        for Positive and ``-1`` for Negative.
        """
        return self.args[2]

    @property
    def sensitivity(self):
        r"""
        Returns the sensitivity function matrix of the feedback loop.

        Sensitivity of a closed-loop system is the ratio of change
        in the open loop gain to the change in the closed loop gain.

        .. note::
            This method would not return the complementary
            sensitivity function.

        Examples
        ========

        >>> from sympy import pprint
        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, MIMOFeedback
        >>> tf1 = TransferFunction(p, 1 - p, p)
        >>> tf2 = TransferFunction(1, p, p)
        >>> tf3 = TransferFunction(1, 1, p)
        >>> sys1 = TransferFunctionMatrix([[tf1, tf2], [tf2, tf1]])
        >>> sys2 = TransferFunctionMatrix([[tf1, tf3], [tf3, tf2]])
        >>> F_1 = MIMOFeedback(sys1, sys2, 1)  # Positive feedback
        >>> F_2 = MIMOFeedback(sys1, sys2)  # Negative feedback
        >>> pprint(F_1.sensitivity, use_unicode=False)
        [   4      3      2               5      4      2           ]
        [- p  + 3*p  - 4*p  + 3*p - 1    p  - 2*p  + 3*p  - 3*p + 1 ]
        [----------------------------  -----------------------------]
        [  4      3      2              5      4      3      2      ]
        [ p  + 3*p  - 8*p  + 8*p - 3   p  + 3*p  - 8*p  + 8*p  - 3*p]
        [                                                           ]
        [       4    3    2                  3      2               ]
        [      p  - p  - p  + p           3*p  - 6*p  + 4*p - 1     ]
        [ --------------------------    --------------------------  ]
        [  4      3      2               4      3      2            ]
        [ p  + 3*p  - 8*p  + 8*p - 3    p  + 3*p  - 8*p  + 8*p - 3  ]
        >>> pprint(F_2.sensitivity, use_unicode=False)
        [ 4      3      2           5      4      2          ]
        [p  - 3*p  + 2*p  + p - 1  p  - 2*p  + 3*p  - 3*p + 1]
        [------------------------  --------------------------]
        [   4      3                   5      4      2       ]
        [  p  - 3*p  + 2*p - 1        p  - 3*p  + 2*p  - p   ]
        [                                                    ]
        [     4    3    2               4      3             ]
        [    p  - p  - p  + p        2*p  - 3*p  + 2*p - 1   ]
        [  -------------------       ---------------------   ]
        [   4      3                   4      3              ]
        [  p  - 3*p  + 2*p - 1        p  - 3*p  + 2*p - 1    ]

        """
        _sys1_mat = self.sys1.doit()._expr_mat
        _sys2_mat = self.sys2.doit()._expr_mat

        return (eye(self.sys1.num_inputs) - \
            self.sign*_sys1_mat*_sys2_mat).inv()

    def doit(self, cancel=True, expand=False, **hints):
        r"""
        Returns the resultant transfer function matrix obtained by the
        feedback interconnection.

        Examples
        ========

        >>> from sympy import pprint
        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, MIMOFeedback
        >>> tf1 = TransferFunction(s, 1 - s, s)
        >>> tf2 = TransferFunction(1, s, s)
        >>> tf3 = TransferFunction(5, 1, s)
        >>> tf4 = TransferFunction(s - 1, s, s)
        >>> tf5 = TransferFunction(0, 1, s)
        >>> sys1 = TransferFunctionMatrix([[tf1, tf2], [tf3, tf4]])
        >>> sys2 = TransferFunctionMatrix([[tf3, tf5], [tf5, tf5]])
        >>> F_1 = MIMOFeedback(sys1, sys2, 1)
        >>> pprint(F_1, use_unicode=False)
        /    [  s      1  ]    [5  0]   \-1   [  s      1  ]
        |    [-----    -  ]    [-  -]   |     [-----    -  ]
        |    [1 - s    s  ]    [1  1]   |     [1 - s    s  ]
        |I - [            ]   *[    ]   |   * [            ]
        |    [  5    s - 1]    [0  0]   |     [  5    s - 1]
        |    [  -    -----]    [-  -]   |     [  -    -----]
        \    [  1      s  ]{t} [1  1]{t}/     [  1      s  ]{t}
        >>> pprint(F_1.doit(), use_unicode=False)
        [  -s           s - 1       ]
        [-------     -----------    ]
        [6*s - 1     s*(6*s - 1)    ]
        [                           ]
        [5*s - 5  (s - 1)*(6*s + 24)]
        [-------  ------------------]
        [6*s - 1     s*(6*s - 1)    ]{t}

        If the user wants the resultant ``TransferFunctionMatrix`` object without
        canceling the common factors then the ``cancel`` kwarg should be passed ``False``.

        >>> pprint(F_1.doit(cancel=False), use_unicode=False)
        [           25*s*(1 - s)                          25 - 25*s              ]
        [       --------------------                    --------------           ]
        [       25*(1 - 6*s)*(1 - s)                    25*s*(1 - 6*s)           ]
        [                                                                        ]
        [s*(25*s - 25) + 5*(1 - s)*(6*s - 1)  s*(s - 1)*(6*s - 1) + s*(25*s - 25)]
        [-----------------------------------  -----------------------------------]
        [         (1 - s)*(6*s - 1)                        2                     ]
        [                                                 s *(6*s - 1)           ]{t}

        If the user wants the expanded form of the resultant transfer function matrix,
        the ``expand`` kwarg should be passed as ``True``.

        >>> pprint(F_1.doit(expand=True), use_unicode=False)
        [  -s          s - 1      ]
        [-------      --------    ]
        [6*s - 1         2        ]
        [             6*s  - s    ]
        [                         ]
        [            2            ]
        [5*s - 5  6*s  + 18*s - 24]
        [-------  ----------------]
        [6*s - 1         2        ]
        [             6*s  - s    ]{t}

        """
        _mat = self.sensitivity * self.sys1.doit()._expr_mat

        _resultant_tfm = _to_TFM(_mat, self.var)

        if cancel:
            _resultant_tfm = _resultant_tfm.simplify()

        if expand:
            _resultant_tfm = _resultant_tfm.expand()

        return _resultant_tfm

    def _eval_rewrite_as_TransferFunctionMatrix(self, sys1, sys2, sign, **kwargs):
        return self.doit()

    def __neg__(self):
        return MIMOFeedback(-self.sys1, -self.sys2, self.sign)


def _to_TFM(mat, var):
    """Private method to convert ImmutableMatrix to TransferFunctionMatrix efficiently"""
    to_tf = lambda expr: TransferFunction.from_rational_expression(expr, var)
    arg = [[to_tf(expr) for expr in row] for row in mat.tolist()]
    return TransferFunctionMatrix(arg)


class TransferFunctionMatrix(MIMOLinearTimeInvariant):
    r"""
    A class for representing the MIMO (multiple-input and multiple-output)
    generalization of the SISO (single-input and single-output) transfer function.

    It is a matrix of transfer functions (``TransferFunction``, SISO-``Series`` or SISO-``Parallel``).
    There is only one argument, ``arg`` which is also the compulsory argument.
    ``arg`` is expected to be strictly of the type list of lists
    which holds the transfer functions or reducible to transfer functions.

    Parameters
    ==========

    arg : Nested ``List`` (strictly).
        Users are expected to input a nested list of ``TransferFunction``, ``Series``
        and/or ``Parallel`` objects.

    Examples
    ========

    .. note::
        ``pprint()`` can be used for better visualization of ``TransferFunctionMatrix`` objects.

    >>> from sympy.abc import s, p, a
    >>> from sympy import pprint
    >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, Series, Parallel
    >>> tf_1 = TransferFunction(s + a, s**2 + s + 1, s)
    >>> tf_2 = TransferFunction(p**4 - 3*p + 2, s + p, s)
    >>> tf_3 = TransferFunction(3, s + 2, s)
    >>> tf_4 = TransferFunction(-a + p, 9*s - 9, s)
    >>> tfm_1 = TransferFunctionMatrix([[tf_1], [tf_2], [tf_3]])
    >>> tfm_1
    TransferFunctionMatrix(((TransferFunction(a + s, s**2 + s + 1, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(3, s + 2, s),)))
    >>> tfm_1.var
    s
    >>> tfm_1.num_inputs
    1
    >>> tfm_1.num_outputs
    3
    >>> tfm_1.shape
    (3, 1)
    >>> tfm_1.args
    (((TransferFunction(a + s, s**2 + s + 1, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(3, s + 2, s),)),)
    >>> tfm_2 = TransferFunctionMatrix([[tf_1, -tf_3], [tf_2, -tf_1], [tf_3, -tf_2]])
    >>> tfm_2
    TransferFunctionMatrix(((TransferFunction(a + s, s**2 + s + 1, s), TransferFunction(-3, s + 2, s)), (TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(-a - s, s**2 + s + 1, s)), (TransferFunction(3, s + 2, s), TransferFunction(-p**4 + 3*p - 2, p + s, s))))
    >>> pprint(tfm_2, use_unicode=False)  # pretty-printing for better visualization
    [   a + s           -3       ]
    [ ----------       -----     ]
    [  2               s + 2     ]
    [ s  + s + 1                 ]
    [                            ]
    [ 4                          ]
    [p  - 3*p + 2      -a - s    ]
    [------------    ----------  ]
    [   p + s         2          ]
    [                s  + s + 1  ]
    [                            ]
    [                 4          ]
    [     3        - p  + 3*p - 2]
    [   -----      --------------]
    [   s + 2          p + s     ]{t}

    TransferFunctionMatrix can be transposed, if user wants to switch the input and output transfer functions

    >>> tfm_2.transpose()
    TransferFunctionMatrix(((TransferFunction(a + s, s**2 + s + 1, s), TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(3, s + 2, s)), (TransferFunction(-3, s + 2, s), TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(-p**4 + 3*p - 2, p + s, s))))
    >>> pprint(_, use_unicode=False)
    [             4                          ]
    [  a + s     p  - 3*p + 2        3       ]
    [----------  ------------      -----     ]
    [ 2             p + s          s + 2     ]
    [s  + s + 1                              ]
    [                                        ]
    [                             4          ]
    [   -3          -a - s     - p  + 3*p - 2]
    [  -----      ----------   --------------]
    [  s + 2       2               p + s     ]
    [             s  + s + 1                 ]{t}

    >>> tf_5 = TransferFunction(5, s, s)
    >>> tf_6 = TransferFunction(5*s, (2 + s**2), s)
    >>> tf_7 = TransferFunction(5, (s*(2 + s**2)), s)
    >>> tf_8 = TransferFunction(5, 1, s)
    >>> tfm_3 = TransferFunctionMatrix([[tf_5, tf_6], [tf_7, tf_8]])
    >>> tfm_3
    TransferFunctionMatrix(((TransferFunction(5, s, s), TransferFunction(5*s, s**2 + 2, s)), (TransferFunction(5, s*(s**2 + 2), s), TransferFunction(5, 1, s))))
    >>> pprint(tfm_3, use_unicode=False)
    [    5        5*s  ]
    [    -       ------]
    [    s        2    ]
    [            s  + 2]
    [                  ]
    [    5         5   ]
    [----------    -   ]
    [  / 2    \    1   ]
    [s*\s  + 2/        ]{t}
    >>> tfm_3.var
    s
    >>> tfm_3.shape
    (2, 2)
    >>> tfm_3.num_outputs
    2
    >>> tfm_3.num_inputs
    2
    >>> tfm_3.args
    (((TransferFunction(5, s, s), TransferFunction(5*s, s**2 + 2, s)), (TransferFunction(5, s*(s**2 + 2), s), TransferFunction(5, 1, s))),)

    To access the ``TransferFunction`` at any index in the ``TransferFunctionMatrix``, use the index notation.

    >>> tfm_3[1, 0]  # gives the TransferFunction present at 2nd Row and 1st Col. Similar to that in Matrix classes
    TransferFunction(5, s*(s**2 + 2), s)
    >>> tfm_3[0, 0]  # gives the TransferFunction present at 1st Row and 1st Col.
    TransferFunction(5, s, s)
    >>> tfm_3[:, 0]  # gives the first column
    TransferFunctionMatrix(((TransferFunction(5, s, s),), (TransferFunction(5, s*(s**2 + 2), s),)))
    >>> pprint(_, use_unicode=False)
    [    5     ]
    [    -     ]
    [    s     ]
    [          ]
    [    5     ]
    [----------]
    [  / 2    \]
    [s*\s  + 2/]{t}
    >>> tfm_3[0, :]  # gives the first row
    TransferFunctionMatrix(((TransferFunction(5, s, s), TransferFunction(5*s, s**2 + 2, s)),))
    >>> pprint(_, use_unicode=False)
    [5   5*s  ]
    [-  ------]
    [s   2    ]
    [   s  + 2]{t}

    To negate a transfer function matrix, ``-`` operator can be prepended:

    >>> tfm_4 = TransferFunctionMatrix([[tf_2], [-tf_1], [tf_3]])
    >>> -tfm_4
    TransferFunctionMatrix(((TransferFunction(-p**4 + 3*p - 2, p + s, s),), (TransferFunction(a + s, s**2 + s + 1, s),), (TransferFunction(-3, s + 2, s),)))
    >>> tfm_5 = TransferFunctionMatrix([[tf_1, tf_2], [tf_3, -tf_1]])
    >>> -tfm_5
    TransferFunctionMatrix(((TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(-p**4 + 3*p - 2, p + s, s)), (TransferFunction(-3, s + 2, s), TransferFunction(a + s, s**2 + s + 1, s))))

    ``subs()`` returns the ``TransferFunctionMatrix`` object with the value substituted in the expression. This will not
    mutate your original ``TransferFunctionMatrix``.

    >>> tfm_2.subs(p, 2)  #  substituting p everywhere in tfm_2 with 2.
    TransferFunctionMatrix(((TransferFunction(a + s, s**2 + s + 1, s), TransferFunction(-3, s + 2, s)), (TransferFunction(12, s + 2, s), TransferFunction(-a - s, s**2 + s + 1, s)), (TransferFunction(3, s + 2, s), TransferFunction(-12, s + 2, s))))
    >>> pprint(_, use_unicode=False)
    [  a + s        -3     ]
    [----------    -----   ]
    [ 2            s + 2   ]
    [s  + s + 1            ]
    [                      ]
    [    12        -a - s  ]
    [  -----     ----------]
    [  s + 2      2        ]
    [            s  + s + 1]
    [                      ]
    [    3          -12    ]
    [  -----       -----   ]
    [  s + 2       s + 2   ]{t}
    >>> pprint(tfm_2, use_unicode=False) # State of tfm_2 is unchanged after substitution
    [   a + s           -3       ]
    [ ----------       -----     ]
    [  2               s + 2     ]
    [ s  + s + 1                 ]
    [                            ]
    [ 4                          ]
    [p  - 3*p + 2      -a - s    ]
    [------------    ----------  ]
    [   p + s         2          ]
    [                s  + s + 1  ]
    [                            ]
    [                 4          ]
    [     3        - p  + 3*p - 2]
    [   -----      --------------]
    [   s + 2          p + s     ]{t}

    ``subs()`` also supports multiple substitutions.

    >>> tfm_2.subs({p: 2, a: 1})  # substituting p with 2 and a with 1
    TransferFunctionMatrix(((TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(-3, s + 2, s)), (TransferFunction(12, s + 2, s), TransferFunction(-s - 1, s**2 + s + 1, s)), (TransferFunction(3, s + 2, s), TransferFunction(-12, s + 2, s))))
    >>> pprint(_, use_unicode=False)
    [  s + 1        -3     ]
    [----------    -----   ]
    [ 2            s + 2   ]
    [s  + s + 1            ]
    [                      ]
    [    12        -s - 1  ]
    [  -----     ----------]
    [  s + 2      2        ]
    [            s  + s + 1]
    [                      ]
    [    3          -12    ]
    [  -----       -----   ]
    [  s + 2       s + 2   ]{t}

    Users can reduce the ``Series`` and ``Parallel`` elements of the matrix to ``TransferFunction`` by using
    ``doit()``.

    >>> tfm_6 = TransferFunctionMatrix([[Series(tf_3, tf_4), Parallel(tf_3, tf_4)]])
    >>> tfm_6
    TransferFunctionMatrix(((Series(TransferFunction(3, s + 2, s), TransferFunction(-a + p, 9*s - 9, s)), Parallel(TransferFunction(3, s + 2, s), TransferFunction(-a + p, 9*s - 9, s))),))
    >>> pprint(tfm_6, use_unicode=False)
    [ -a + p   3     -a + p     3  ]
    [-------*-----  ------- + -----]
    [9*s - 9 s + 2  9*s - 9   s + 2]{t}
    >>> tfm_6.doit()
    TransferFunctionMatrix(((TransferFunction(-3*a + 3*p, (s + 2)*(9*s - 9), s), TransferFunction(27*s + (-a + p)*(s + 2) - 27, (s + 2)*(9*s - 9), s)),))
    >>> pprint(_, use_unicode=False)
    [    -3*a + 3*p     27*s + (-a + p)*(s + 2) - 27]
    [-----------------  ----------------------------]
    [(s + 2)*(9*s - 9)       (s + 2)*(9*s - 9)      ]{t}
    >>> tf_9 = TransferFunction(1, s, s)
    >>> tf_10 = TransferFunction(1, s**2, s)
    >>> tfm_7 = TransferFunctionMatrix([[Series(tf_9, tf_10), tf_9], [tf_10, Parallel(tf_9, tf_10)]])
    >>> tfm_7
    TransferFunctionMatrix(((Series(TransferFunction(1, s, s), TransferFunction(1, s**2, s)), TransferFunction(1, s, s)), (TransferFunction(1, s**2, s), Parallel(TransferFunction(1, s, s), TransferFunction(1, s**2, s)))))
    >>> pprint(tfm_7, use_unicode=False)
    [ 1      1   ]
    [----    -   ]
    [   2    s   ]
    [s*s         ]
    [            ]
    [ 1    1    1]
    [ --   -- + -]
    [  2    2   s]
    [ s    s     ]{t}
    >>> tfm_7.doit()
    TransferFunctionMatrix(((TransferFunction(1, s**3, s), TransferFunction(1, s, s)), (TransferFunction(1, s**2, s), TransferFunction(s**2 + s, s**3, s))))
    >>> pprint(_, use_unicode=False)
    [1     1   ]
    [--    -   ]
    [ 3    s   ]
    [s         ]
    [          ]
    [     2    ]
    [1   s  + s]
    [--  ------]
    [ 2     3  ]
    [s     s   ]{t}

    Addition, subtraction, and multiplication of transfer function matrices can form
    unevaluated ``Series`` or ``Parallel`` objects.

    - For addition and subtraction:
      All the transfer function matrices must have the same shape.

    - For multiplication (C = A * B):
      The number of inputs of the first transfer function matrix (A) must be equal to the
      number of outputs of the second transfer function matrix (B).

    Also, use pretty-printing (``pprint``) to analyse better.

    >>> tfm_8 = TransferFunctionMatrix([[tf_3], [tf_2], [-tf_1]])
    >>> tfm_9 = TransferFunctionMatrix([[-tf_3]])
    >>> tfm_10 = TransferFunctionMatrix([[tf_1], [tf_2], [tf_4]])
    >>> tfm_11 = TransferFunctionMatrix([[tf_4], [-tf_1]])
    >>> tfm_12 = TransferFunctionMatrix([[tf_4, -tf_1, tf_3], [-tf_2, -tf_4, -tf_3]])
    >>> tfm_8 + tfm_10
    MIMOParallel(TransferFunctionMatrix(((TransferFunction(3, s + 2, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(-a - s, s**2 + s + 1, s),))), TransferFunctionMatrix(((TransferFunction(a + s, s**2 + s + 1, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(-a + p, 9*s - 9, s),))))
    >>> pprint(_, use_unicode=False)
    [     3      ]      [   a + s    ]
    [   -----    ]      [ ---------- ]
    [   s + 2    ]      [  2         ]
    [            ]      [ s  + s + 1 ]
    [ 4          ]      [            ]
    [p  - 3*p + 2]      [ 4          ]
    [------------]    + [p  - 3*p + 2]
    [   p + s    ]      [------------]
    [            ]      [   p + s    ]
    [   -a - s   ]      [            ]
    [ ---------- ]      [   -a + p   ]
    [  2         ]      [  -------   ]
    [ s  + s + 1 ]{t}   [  9*s - 9   ]{t}
    >>> -tfm_10 - tfm_8
    MIMOParallel(TransferFunctionMatrix(((TransferFunction(-a - s, s**2 + s + 1, s),), (TransferFunction(-p**4 + 3*p - 2, p + s, s),), (TransferFunction(a - p, 9*s - 9, s),))), TransferFunctionMatrix(((TransferFunction(-3, s + 2, s),), (TransferFunction(-p**4 + 3*p - 2, p + s, s),), (TransferFunction(a + s, s**2 + s + 1, s),))))
    >>> pprint(_, use_unicode=False)
    [    -a - s    ]      [     -3       ]
    [  ----------  ]      [    -----     ]
    [   2          ]      [    s + 2     ]
    [  s  + s + 1  ]      [              ]
    [              ]      [   4          ]
    [   4          ]      [- p  + 3*p - 2]
    [- p  + 3*p - 2]    + [--------------]
    [--------------]      [    p + s     ]
    [    p + s     ]      [              ]
    [              ]      [    a + s     ]
    [    a - p     ]      [  ----------  ]
    [   -------    ]      [   2          ]
    [   9*s - 9    ]{t}   [  s  + s + 1  ]{t}
    >>> tfm_12 * tfm_8
    MIMOSeries(TransferFunctionMatrix(((TransferFunction(3, s + 2, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(-a - s, s**2 + s + 1, s),))), TransferFunctionMatrix(((TransferFunction(-a + p, 9*s - 9, s), TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(3, s + 2, s)), (TransferFunction(-p**4 + 3*p - 2, p + s, s), TransferFunction(a - p, 9*s - 9, s), TransferFunction(-3, s + 2, s)))))
    >>> pprint(_, use_unicode=False)
                                           [     3      ]
                                           [   -----    ]
    [    -a + p        -a - s      3  ]    [   s + 2    ]
    [   -------      ----------  -----]    [            ]
    [   9*s - 9       2          s + 2]    [ 4          ]
    [                s  + s + 1       ]    [p  - 3*p + 2]
    [                                 ]   *[------------]
    [   4                             ]    [   p + s    ]
    [- p  + 3*p - 2    a - p      -3  ]    [            ]
    [--------------   -------    -----]    [   -a - s   ]
    [    p + s        9*s - 9    s + 2]{t} [ ---------- ]
                                           [  2         ]
                                           [ s  + s + 1 ]{t}
    >>> tfm_12 * tfm_8 * tfm_9
    MIMOSeries(TransferFunctionMatrix(((TransferFunction(-3, s + 2, s),),)), TransferFunctionMatrix(((TransferFunction(3, s + 2, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(-a - s, s**2 + s + 1, s),))), TransferFunctionMatrix(((TransferFunction(-a + p, 9*s - 9, s), TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(3, s + 2, s)), (TransferFunction(-p**4 + 3*p - 2, p + s, s), TransferFunction(a - p, 9*s - 9, s), TransferFunction(-3, s + 2, s)))))
    >>> pprint(_, use_unicode=False)
                                           [     3      ]
                                           [   -----    ]
    [    -a + p        -a - s      3  ]    [   s + 2    ]
    [   -------      ----------  -----]    [            ]
    [   9*s - 9       2          s + 2]    [ 4          ]
    [                s  + s + 1       ]    [p  - 3*p + 2]    [ -3  ]
    [                                 ]   *[------------]   *[-----]
    [   4                             ]    [   p + s    ]    [s + 2]{t}
    [- p  + 3*p - 2    a - p      -3  ]    [            ]
    [--------------   -------    -----]    [   -a - s   ]
    [    p + s        9*s - 9    s + 2]{t} [ ---------- ]
                                           [  2         ]
                                           [ s  + s + 1 ]{t}
    >>> tfm_10 + tfm_8*tfm_9
    MIMOParallel(TransferFunctionMatrix(((TransferFunction(a + s, s**2 + s + 1, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(-a + p, 9*s - 9, s),))), MIMOSeries(TransferFunctionMatrix(((TransferFunction(-3, s + 2, s),),)), TransferFunctionMatrix(((TransferFunction(3, s + 2, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(-a - s, s**2 + s + 1, s),)))))
    >>> pprint(_, use_unicode=False)
    [   a + s    ]      [     3      ]
    [ ---------- ]      [   -----    ]
    [  2         ]      [   s + 2    ]
    [ s  + s + 1 ]      [            ]
    [            ]      [ 4          ]
    [ 4          ]      [p  - 3*p + 2]    [ -3  ]
    [p  - 3*p + 2]    + [------------]   *[-----]
    [------------]      [   p + s    ]    [s + 2]{t}
    [   p + s    ]      [            ]
    [            ]      [   -a - s   ]
    [   -a + p   ]      [ ---------- ]
    [  -------   ]      [  2         ]
    [  9*s - 9   ]{t}   [ s  + s + 1 ]{t}

    These unevaluated ``Series`` or ``Parallel`` objects can convert into the
    resultant transfer function matrix using ``.doit()`` method or by
    ``.rewrite(TransferFunctionMatrix)``.

    >>> (-tfm_8 + tfm_10 + tfm_8*tfm_9).doit()
    TransferFunctionMatrix(((TransferFunction((a + s)*(s + 2)**3 - 3*(s + 2)**2*(s**2 + s + 1) - 9*(s + 2)*(s**2 + s + 1), (s + 2)**3*(s**2 + s + 1), s),), (TransferFunction((p + s)*(-3*p**4 + 9*p - 6), (p + s)**2*(s + 2), s),), (TransferFunction((-a + p)*(s + 2)*(s**2 + s + 1)**2 + (a + s)*(s + 2)*(9*s - 9)*(s**2 + s + 1) + (3*a + 3*s)*(9*s - 9)*(s**2 + s + 1), (s + 2)*(9*s - 9)*(s**2 + s + 1)**2, s),)))
    >>> (-tfm_12 * -tfm_8 * -tfm_9).rewrite(TransferFunctionMatrix)
    TransferFunctionMatrix(((TransferFunction(3*(-3*a + 3*p)*(p + s)*(s + 2)*(s**2 + s + 1)**2 + 3*(-3*a - 3*s)*(p + s)*(s + 2)*(9*s - 9)*(s**2 + s + 1) + 3*(a + s)*(s + 2)**2*(9*s - 9)*(-p**4 + 3*p - 2)*(s**2 + s + 1), (p + s)*(s + 2)**3*(9*s - 9)*(s**2 + s + 1)**2, s),), (TransferFunction(3*(-a + p)*(p + s)*(s + 2)**2*(-p**4 + 3*p - 2)*(s**2 + s + 1) + 3*(3*a + 3*s)*(p + s)**2*(s + 2)*(9*s - 9) + 3*(p + s)*(s + 2)*(9*s - 9)*(-3*p**4 + 9*p - 6)*(s**2 + s + 1), (p + s)**2*(s + 2)**3*(9*s - 9)*(s**2 + s + 1), s),)))

    See Also
    ========

    TransferFunction, MIMOSeries, MIMOParallel, Feedback

    """
    def __new__(cls, arg):

        expr_mat_arg = []
        try:
            var = arg[0][0].var
        except TypeError:
            raise ValueError("`arg` param in TransferFunctionMatrix should "
            "strictly be a nested list containing TransferFunction objects.")
        for row_index, row in enumerate(arg):
            temp = []
            for col_index, element in enumerate(row):
                if not isinstance(element, SISOLinearTimeInvariant):
                    raise TypeError("Each element is expected to be of type `SISOLinearTimeInvariant`.")

                if var != element.var:
                    raise ValueError("Conflicting value(s) found for `var`. All TransferFunction instances in "
                            "TransferFunctionMatrix should use the same complex variable in Laplace domain.")

                temp.append(element.to_expr())
            expr_mat_arg.append(temp)

        if isinstance(arg, (tuple, list, Tuple)):
            # Making nested Tuple (sympy.core.containers.Tuple) from nested list or nested Python tuple
            arg = Tuple(*(Tuple(*r, sympify=False) for r in arg), sympify=False)

        obj = super(TransferFunctionMatrix, cls).__new__(cls, arg)
        obj._expr_mat = ImmutableMatrix(expr_mat_arg)
        return obj

    @classmethod
    def from_Matrix(cls, matrix, var):
        """
        Creates a new ``TransferFunctionMatrix`` efficiently from a SymPy Matrix of ``Expr`` objects.

        Parameters
        ==========

        matrix : ``ImmutableMatrix`` having ``Expr``/``Number`` elements.
        var : Symbol
            Complex variable of the Laplace transform which will be used by the
            all the ``TransferFunction`` objects in the ``TransferFunctionMatrix``.

        Examples
        ========

        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunctionMatrix
        >>> from sympy import Matrix, pprint
        >>> M = Matrix([[s, 1/s], [1/(s+1), s]])
        >>> M_tf = TransferFunctionMatrix.from_Matrix(M, s)
        >>> pprint(M_tf, use_unicode=False)
        [  s    1]
        [  -    -]
        [  1    s]
        [        ]
        [  1    s]
        [-----  -]
        [s + 1  1]{t}
        >>> M_tf.elem_poles()
        [[[], [0]], [[-1], []]]
        >>> M_tf.elem_zeros()
        [[[0], []], [[], [0]]]

        """
        return _to_TFM(matrix, var)

    @property
    def var(self):
        """
        Returns the complex variable used by all the transfer functions or
        ``Series``/``Parallel`` objects in a transfer function matrix.

        Examples
        ========

        >>> from sympy.abc import p, s
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, Series, Parallel
        >>> G1 = TransferFunction(p**2 + 2*p + 4, p - 6, p)
        >>> G2 = TransferFunction(p, 4 - p, p)
        >>> G3 = TransferFunction(0, p**4 - 1, p)
        >>> G4 = TransferFunction(s + 1, s**2 + s + 1, s)
        >>> S1 = Series(G1, G2)
        >>> S2 = Series(-G3, Parallel(G2, -G1))
        >>> tfm1 = TransferFunctionMatrix([[G1], [G2], [G3]])
        >>> tfm1.var
        p
        >>> tfm2 = TransferFunctionMatrix([[-S1, -S2], [S1, S2]])
        >>> tfm2.var
        p
        >>> tfm3 = TransferFunctionMatrix([[G4]])
        >>> tfm3.var
        s

        """
        return self.args[0][0][0].var

    @property
    def num_inputs(self):
        """
        Returns the number of inputs of the system.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
        >>> G1 = TransferFunction(s + 3, s**2 - 3, s)
        >>> G2 = TransferFunction(4, s**2, s)
        >>> G3 = TransferFunction(p**2 + s**2, p - 3, s)
        >>> tfm_1 = TransferFunctionMatrix([[G2, -G1, G3], [-G2, -G1, -G3]])
        >>> tfm_1.num_inputs
        3

        See Also
        ========

        num_outputs

        """
        return self._expr_mat.shape[1]

    @property
    def num_outputs(self):
        """
        Returns the number of outputs of the system.

        Examples
        ========

        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunctionMatrix
        >>> from sympy import Matrix
        >>> M_1 = Matrix([[s], [1/s]])
        >>> TFM = TransferFunctionMatrix.from_Matrix(M_1, s)
        >>> print(TFM)
        TransferFunctionMatrix(((TransferFunction(s, 1, s),), (TransferFunction(1, s, s),)))
        >>> TFM.num_outputs
        2

        See Also
        ========

        num_inputs

        """
        return self._expr_mat.shape[0]

    @property
    def shape(self):
        """
        Returns the shape of the transfer function matrix, that is, ``(# of outputs, # of inputs)``.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
        >>> tf1 = TransferFunction(p**2 - 1, s**4 + s**3 - p, p)
        >>> tf2 = TransferFunction(1 - p, p**2 - 3*p + 7, p)
        >>> tf3 = TransferFunction(3, 4, p)
        >>> tfm1 = TransferFunctionMatrix([[tf1, -tf2]])
        >>> tfm1.shape
        (1, 2)
        >>> tfm2 = TransferFunctionMatrix([[-tf2, tf3], [tf1, -tf1]])
        >>> tfm2.shape
        (2, 2)

        """
        return self._expr_mat.shape

    def __neg__(self):
        neg = -self._expr_mat
        return _to_TFM(neg, self.var)

    @_check_other_MIMO
    def __add__(self, other):

        if not isinstance(other, MIMOParallel):
            return MIMOParallel(self, other)
        other_arg_list = list(other.args)
        return MIMOParallel(self, *other_arg_list)

    @_check_other_MIMO
    def __sub__(self, other):
        return self + (-other)

    @_check_other_MIMO
    def __mul__(self, other):

        if not isinstance(other, MIMOSeries):
            return MIMOSeries(other, self)
        other_arg_list = list(other.args)
        return MIMOSeries(*other_arg_list, self)

    def __getitem__(self, key):
        trunc = self._expr_mat.__getitem__(key)
        if isinstance(trunc, ImmutableMatrix):
            return _to_TFM(trunc, self.var)
        return TransferFunction.from_rational_expression(trunc, self.var)

    def transpose(self):
        """Returns the transpose of the ``TransferFunctionMatrix`` (switched input and output layers)."""
        transposed_mat = self._expr_mat.transpose()
        return _to_TFM(transposed_mat, self.var)

    def elem_poles(self):
        """
        Returns the poles of each element of the ``TransferFunctionMatrix``.

        .. note::
            Actual poles of a MIMO system are NOT the poles of individual elements.

        Examples
        ========

        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
        >>> tf_1 = TransferFunction(3, (s + 1), s)
        >>> tf_2 = TransferFunction(s + 6, (s + 1)*(s + 2), s)
        >>> tf_3 = TransferFunction(s + 3, s**2 + 3*s + 2, s)
        >>> tf_4 = TransferFunction(s + 2, s**2 + 5*s - 10, s)
        >>> tfm_1 = TransferFunctionMatrix([[tf_1, tf_2], [tf_3, tf_4]])
        >>> tfm_1
        TransferFunctionMatrix(((TransferFunction(3, s + 1, s), TransferFunction(s + 6, (s + 1)*(s + 2), s)), (TransferFunction(s + 3, s**2 + 3*s + 2, s), TransferFunction(s + 2, s**2 + 5*s - 10, s))))
        >>> tfm_1.elem_poles()
        [[[-1], [-2, -1]], [[-2, -1], [-5/2 + sqrt(65)/2, -sqrt(65)/2 - 5/2]]]

        See Also
        ========

        elem_zeros

        """
        return [[element.poles() for element in row] for row in self.doit().args[0]]

    def elem_zeros(self):
        """
        Returns the zeros of each element of the ``TransferFunctionMatrix``.

        .. note::
            Actual zeros of a MIMO system are NOT the zeros of individual elements.

        Examples
        ========

        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
        >>> tf_1 = TransferFunction(3, (s + 1), s)
        >>> tf_2 = TransferFunction(s + 6, (s + 1)*(s + 2), s)
        >>> tf_3 = TransferFunction(s + 3, s**2 + 3*s + 2, s)
        >>> tf_4 = TransferFunction(s**2 - 9*s + 20, s**2 + 5*s - 10, s)
        >>> tfm_1 = TransferFunctionMatrix([[tf_1, tf_2], [tf_3, tf_4]])
        >>> tfm_1
        TransferFunctionMatrix(((TransferFunction(3, s + 1, s), TransferFunction(s + 6, (s + 1)*(s + 2), s)), (TransferFunction(s + 3, s**2 + 3*s + 2, s), TransferFunction(s**2 - 9*s + 20, s**2 + 5*s - 10, s))))
        >>> tfm_1.elem_zeros()
        [[[], [-6]], [[-3], [4, 5]]]

        See Also
        ========

        elem_poles

        """
        return [[element.zeros() for element in row] for row in self.doit().args[0]]

    def _flat(self):
        """Returns flattened list of args in TransferFunctionMatrix"""
        return [elem for tup in self.args[0] for elem in tup]

    def _eval_evalf(self, prec):
        """Calls evalf() on each transfer function in the transfer function matrix"""
        dps = prec_to_dps(prec)
        mat = self._expr_mat.applyfunc(lambda a: a.evalf(n=dps))
        return _to_TFM(mat, self.var)

    def _eval_simplify(self, **kwargs):
        """Simplifies the transfer function matrix"""
        simp_mat = self._expr_mat.applyfunc(lambda a: cancel(a, expand=False))
        return _to_TFM(simp_mat, self.var)

    def expand(self, **hints):
        """Expands the transfer function matrix"""
        expand_mat = self._expr_mat.expand(**hints)
        return _to_TFM(expand_mat, self.var)
