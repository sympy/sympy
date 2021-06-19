from sympy import Basic, Mul, Pow, degree, Symbol, expand, cancel, Expr, exp, roots
from sympy.core.evalf import EvalfMixin
from sympy.core.logic import fuzzy_and
from sympy.core.numbers import Integer
from sympy.core.sympify import sympify, _sympify
from sympy.polys import Poly, rootof
from sympy.series import limit

__all__ = ['TransferFunction', 'Series', 'Parallel', 'Feedback']


def _roots(poly, var):
    """ like roots, but works on higher-order polynomials. """
    r = roots(poly, var, multiple=True)
    n = degree(poly)
    if len(r) != n:
        r = [rootof(poly, var, k) for k in range(n)]
    return r


class TransferFunction(Basic, EvalfMixin):
    r"""
    A class for representing LTI (Linear, time-invariant) systems that can be strictly described
    by ratio of polynomials in the Laplace transform complex variable. The arguments
    are ``num``, ``den``, and ``var``, where ``num`` and ``den`` are numerator and
    denominator polynomials of the ``TransferFunction`` respectively, and the third argument is
    a complex variable of the Laplace transform used by these polynomials of the transfer function.
    ``num`` and ``den`` can be either polynomials or numbers, whereas ``var``
    has to be a Symbol.

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

    The numerator of the transfer function is therefore, the Laplace transform of the output signal
    (The signals are represented as functions of time) and similarly the denominator
    of the transfer function is the Laplace transform of the input signal. It is also a convention
    to denote the input and output signal's Laplace transform with capital alphabets like shown below.

            $H(s) = \frac{Y(s)}{X(s)} = \frac{ \mathcal{L}\left\{y(t)\right\} }{ \mathcal{L}\left\{x(t)\right\} }$

    $s$, also known as complex frequency, is a complex variable in the the Laplace domain. It corresponds to the
    equivalent variable $t$, in the time domain. Independent variable $t$ gets transformed to $s$ after Laplace
    transform. To get back to the time domain from Laplace domain, Inverse-Laplace transform is used.
    Transfer function, $H$, is generally given as a rational function in $s$ as-

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
    ZeroDivisionError
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

    >>> tf2 = TransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p)
    >>> tf2
    TransferFunction(p**2 + 2*p - 3, p**2 + 4*p - 5, p)

    To negate a transfer function the ``-`` operator can be prepended:

    >>> tf3 = TransferFunction(-a + s, p**2 + s, p)
    >>> -tf3
    TransferFunction(a - s, p**2 + s, p)
    >>> tf4 = TransferFunction(s**4 - 2*s**3 + 5*s + 4, s + 4, s)
    >>> -tf4
    TransferFunction(-s**4 + 2*s**3 - 5*s - 4, s + 4, s)

    You can use a Float or an Integer (or other constants) as numerator and denominator:

    >>> tf5 = TransferFunction(1/2, 4, s)
    >>> tf5.num
    0.500000000000000
    >>> tf5.den
    4
    >>> tf5.var
    s
    >>> tf5.args
    (0.5, 4, s)

    You can take the integer power of a transfer function using the ``**`` operator:

    >>> tf6 = TransferFunction(s + a, s - a, s)
    >>> tf6**3
    TransferFunction(a**3 + 3*a**2*s + 3*a*s**2 + s**3, -a**3 + 3*a**2*s - 3*a*s**2 + s**3, s)
    >>> tf6**0
    TransferFunction(1, 1, s)
    >>> tf7 = TransferFunction(p + 4, p - 3, p)
    >>> tf7**-1
    TransferFunction(p - 3, p + 4, p)

    Addition, subtraction, and multiplication of transfer functions can form
    unevaluated ``Series`` or ``Parallel`` objects.

    >>> tf8 = TransferFunction(s + 1, s**2 + s + 1, s)
    >>> tf9 = TransferFunction(s - p, s + 3, s)
    >>> tf10 = TransferFunction(4*s**2 + 2*s - 4, s - 1, s)
    >>> tf11 = TransferFunction(1 - s, s**2 + 4, s)
    >>> tf8 + tf9
    Parallel(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(-p + s, s + 3, s))
    >>> tf9 - tf10
    Parallel(TransferFunction(-p + s, s + 3, s), TransferFunction(-4*s**2 - 2*s + 4, s - 1, s))
    >>> tf8 * tf9
    Series(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(-p + s, s + 3, s))
    >>> tf9 - (tf8 + tf11)
    Parallel(TransferFunction(-p + s, s + 3, s), TransferFunction(-s - 1, s**2 + s + 1, s), TransferFunction(s - 1, s**2 + 4, s))
    >>> tf9 - (tf8 * tf11)
    Parallel(TransferFunction(-p + s, s + 3, s), Series(TransferFunction(-1, 1, s), Series(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(1 - s, s**2 + 4, s))))
    >>> tf10 * tf9 * tf8
    Series(TransferFunction(4*s**2 + 2*s - 4, s - 1, s), TransferFunction(-p + s, s + 3, s), TransferFunction(s + 1, s**2 + s + 1, s))
    >>> tf8 * tf10 + tf9 * tf11
    Parallel(Series(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(4*s**2 + 2*s - 4, s - 1, s)), Series(TransferFunction(-p + s, s + 3, s), TransferFunction(1 - s, s**2 + 4, s)))
    >>> (tf8 + tf11) * (tf9 + tf10)
    Series(Parallel(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(1 - s, s**2 + 4, s)), Parallel(TransferFunction(-p + s, s + 3, s), TransferFunction(4*s**2 + 2*s - 4, s - 1, s)))

    These unevaluated ``Series`` or ``Parallel`` objects can convert into the
    resultant transfer function using ``.doit()`` method or by ``.rewrite(TransferFunction)``.

    >>> ((tf8 + tf9) * tf11).doit()
    TransferFunction(p*s**3 - p - s**4 - s**3 - 3*s**2 + 2*s + 3, s**5 + 4*s**4 + 8*s**3 + 19*s**2 + 16*s + 12, s)
    >>> (tf8 * tf9 - tf10 * tf11).rewrite(TransferFunction)
    TransferFunction(-p*s**4 - 3*p*s**2 + 4*p + 4*s**6 + 15*s**5 + 2*s**4 - 13*s**3 - 14*s**2 - 6*s + 12, s**6 + 3*s**5 + 4*s**4 + 11*s**3 - 3*s**2 - 4*s - 12, s)

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
        num, den = expand(num), expand(den)

        if not isinstance(var, Symbol):
            raise TypeError("Variable input must be a Symbol.")

        if den == 0:
            raise ZeroDivisionError("TransferFunction can't have a zero denominator.")

        if (((isinstance(num, Expr) and num.has(Symbol) and not num.has(exp)) or num.is_number) and
            ((isinstance(den, Expr) and den.has(Symbol) and not den.has(exp)) or den.is_number)):
            obj = super(TransferFunction, cls).__new__(cls, num, den, var)
            obj._num = num
            obj._den = den
            obj._var = var
            return obj

        else:
            raise TypeError("Unsupported type for numerator or denominator of TransferFunction.")

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
        p**2 + 2*p - 15

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
        >>> from sympy.core.symbol import symbols
        >>> q, r = symbols('q, r', negative=True)
        >>> from sympy.physics.control.lti import TransferFunction
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
            raise TypeError("TransferFunction cannot be added with {}.".
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
            raise TypeError("{} cannot be subtracted from a TransferFunction."
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
            raise TypeError("TransferFunction cannot be multiplied with {}."
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
            raise TypeError("TransferFunction cannot be divided by {}.".
                format(type(other)))

    __rtruediv__ = __truediv__

    def __pow__(self, p):
        p = sympify(p)
        if not isinstance(p, Integer):
            raise ValueError("Exponent must be an Integer.")
        if p == 0:
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

    def _to_expr(self):
        """To convert TransferFunction type to SymPy Expr"""
        return Mul(self.num, Pow(self.den, -1, evaluate=False), evaluate=False)


class Series(Basic):
    """
    A class for representing product of transfer functions or transfer functions in a
    series configuration.

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
    TransferFunction(-a*p**4*s**3 + 2*a*p**4 - b*p**2*s**4 + 2*b*p**2*s, -p**2*s**4 - 5*p**2*s - 6*p**2 + s**6 + 5*s**3 + 6*s**2, s)
    >>> S4 = Series(tf2, Parallel(tf1, -tf3))
    >>> S4.doit()
    TransferFunction(a*p**3*s**3 - 2*a*p**3 + a*p**2*s**4 - 2*a*p**2*s + b*p*s**4 - 2*b*p*s + b*s**5 - 2*b*s**2 + p**3*s**3 - 2*p**3 - p**2*s**4 + 2*p**2*s, -p**2*s**4 - 5*p**2*s - 6*p**2 + s**6 + 5*s**3 + 6*s**2, s)

    Notes
    =====

    All the transfer functions should use the same complex variable
    ``var`` of the Laplace transform.

    See Also
    ========

    Parallel, TransferFunction, Feedback

    """
    def __new__(cls, *args, evaluate=False):
        if not all(isinstance(arg, (TransferFunction, Parallel, Series)) for arg in args):
            raise TypeError("Unsupported type of argument(s) for Series.")

        obj = super().__new__(cls, *args)
        obj._var = None
        for arg in args:
            if obj._var is None:
                obj._var = arg.var
            elif obj._var != arg.var:
                raise ValueError("All transfer functions should use the same complex"
                    " variable of the Laplace transform.")
        if evaluate:
            return obj.doit()
        return obj

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
        return self._var

    def doit(self, **kwargs):
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
        TransferFunction(a*p**2*s**3 - 2*a*p**2 + b*s**4 - 2*b*s, -p*s**4 - 5*p*s - 6*p + s**5 + 5*s**2 + 6*s, s)
        >>> Series(-tf1, -tf2).doit()
        TransferFunction(a*p**2*s**3 - 2*a*p**2 + b*s**4 - 2*b*s, -p*s**4 - 5*p*s - 6*p + s**5 + 5*s**2 + 6*s, s)

        """
        res = None
        for arg in self.args:
            arg = arg.doit()
            if res is None:
                res = arg
            else:
                num_ = arg.num * res.num
                den_ = arg.den * res.den
                res = TransferFunction(num_, den_, self.var)
        return res

    def _eval_rewrite_as_TransferFunction(self, *args, **kwargs):
        return self.doit()

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
            raise ValueError("This transfer function expression is invalid.")

    __radd__ = __add__

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
            raise ValueError("This transfer function expression is invalid.")

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if isinstance(other, (TransferFunction, Parallel)):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            arg_list = list(self.args)

            return Series(*arg_list, other)
        elif isinstance(other, Series):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            self_arg_list = list(self.args)
            other_arg_list = list(other.args)

            return Series(*self_arg_list, *other_arg_list)
        else:
            raise ValueError("This transfer function expression is invalid.")

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


class Parallel(Basic):
    """
    A class for representing addition of transfer functions or transfer functions
    in a parallel configuration.

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
    TransferFunction(a*p**3*s**4 + 5*a*p**3*s + 6*a*p**3 + a*p**2*s**5 + 5*a*p**2*s**2 + 6*a*p**2*s + b*p*s**5 + 5*b*p*s**2 + 6*b*p*s + b*s**6 + 5*b*s**3 + 6*b*s**2 + p**3*s**4 + 5*p**3*s + 6*p**3 - p**2*s**5 - p**2*s**3 - 5*p**2*s**2 - 6*p**2*s + 2*p**2 + s**5 - 2*s**2, -p**2*s**4 - 5*p**2*s - 6*p**2 + s**6 + 5*s**3 + 6*s**2, s)
    >>> Parallel(tf2, Series(tf1, -tf3)).doit()
    TransferFunction(-a*p**4*s**4 - 5*a*p**4*s - 6*a*p**4 - b*p**2*s**5 - 5*b*p**2*s**2 - 6*b*p**2*s - p**2*s**3 + 2*p**2 + s**5 - 2*s**2, -p**2*s**4 - 5*p**2*s - 6*p**2 + s**6 + 5*s**3 + 6*s**2, s)

    Notes
    =====

    All the transfer functions should use the same complex variable
    ``var`` of the Laplace transform.

    See Also
    ========

    Series, TransferFunction, Feedback

    """
    def __new__(cls, *args, evaluate=False):
        if not all(isinstance(arg, (TransferFunction, Series, Parallel)) for arg in args):
            raise TypeError("Unsupported type of argument(s) for Parallel.")

        obj = super().__new__(cls, *args)
        obj._var = None
        for arg in args:
            if obj._var is None:
                obj._var = arg.var
            elif obj._var != arg.var:
                raise ValueError("All transfer functions should use the same complex"
                    " variable of the Laplace transform.")
        if evaluate:
            return obj.doit()
        return obj

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
        return self._var

    def doit(self, **kwargs):
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
        TransferFunction(a*p**2*s**4 + 5*a*p**2*s + 6*a*p**2 + b*s**5 + 5*b*s**2 + 6*b*s - p*s**3 + 2*p + s**4 - 2*s, -p*s**4 - 5*p*s - 6*p + s**5 + 5*s**2 + 6*s, s)
        >>> Parallel(-tf1, -tf2).doit()
        TransferFunction(-a*p**2*s**4 - 5*a*p**2*s - 6*a*p**2 - b*s**5 - 5*b*s**2 - 6*b*s + p*s**3 - 2*p - s**4 + 2*s, -p*s**4 - 5*p*s - 6*p + s**5 + 5*s**2 + 6*s, s)

        """
        res = None
        for arg in self.args:
            arg = arg.doit()
            if res is None:
                res = arg
            else:
                num_ = res.num * arg.den + res.den * arg.num
                den_ = res.den * arg.den
                res = TransferFunction(num_, den_, self.var)
        return res

    def _eval_rewrite_as_TransferFunction(self, *args, **kwargs):
        return self.doit()

    def __add__(self, other):
        if isinstance(other, (TransferFunction, Series)):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            arg_list = list(self.args)
            arg_list.append(other)

            return Parallel(*arg_list)
        elif isinstance(other, Parallel):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            self_arg_list = list(self.args)
            other_arg_list = list(other.args)
            for elem in other_arg_list:
                self_arg_list.append(elem)

            return Parallel(*self_arg_list)
        else:
            raise ValueError("This transfer function expression is invalid.")

    def __sub__(self, other):
        if isinstance(other, (TransferFunction, Series)):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            arg_list = list(self.args)
            arg_list.append(-other)

            return Parallel(*arg_list)
        elif isinstance(other, Parallel):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            self_arg_list = list(self.args)
            other_arg_list = list(other.args)
            for elem in other_arg_list:
                self_arg_list.append(-elem)

            return Parallel(*self_arg_list)
        else:
            raise ValueError("This transfer function expression is invalid.")

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
            raise ValueError("This transfer function expression is invalid.")

    def __neg__(self):
        return Series(TransferFunction(-1, 1, self.var), self)

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


class Feedback(Basic):
    """
    A class for representing negative feedback interconnection between two
    input/output systems. The first argument, ``num``, is called as the
    primary plant or the numerator, and the second argument, ``den``, is
    called as the feedback plant (which is often a feedback controller) or
    the denominator. Both ``num`` and ``den`` can either be ``Series`` or
    ``TransferFunction`` objects.

    Parameters
    ==========

    num : Series, TransferFunction
        The primary plant.
    den : Series, TransferFunction
        The feedback plant (often a feedback controller).

    Raises
    ======

    ValueError
        When ``num`` is equal to ``den`` or when they are not using the
        same complex variable of the Laplace transform.
    TypeError
        When either ``num`` or ``den`` is not a ``Series`` or a
        ``TransferFunction`` object.

    Examples
    ========

    >>> from sympy.abc import s
    >>> from sympy.physics.control.lti import TransferFunction, Feedback
    >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
    >>> controller = TransferFunction(5*s - 10, s + 7, s)
    >>> F1 = Feedback(plant, controller)
    >>> F1
    Feedback(TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s), TransferFunction(5*s - 10, s + 7, s))
    >>> F1.var
    s
    >>> F1.args
    (TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s), TransferFunction(5*s - 10, s + 7, s))

    You can get the primary and the feedback plant using ``.num`` and ``.den`` respectively.

    >>> F1.num
    TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
    >>> F1.den
    TransferFunction(5*s - 10, s + 7, s)

    You can get the resultant closed loop transfer function obtained by negative feedback
    interconnection using ``.doit()`` method.

    >>> F1.doit()
    TransferFunction(3*s**5 + 16*s**4 - 60*s**3 - 149*s**2 + 176*s - 42, 16*s**5 - 56*s**4 - 111*s**3 + 504*s**2 - 398*s + 88, s)
    >>> G = TransferFunction(2*s**2 + 5*s + 1, s**2 + 2*s + 3, s)
    >>> C = TransferFunction(5*s + 10, s + 10, s)
    >>> F2 = Feedback(G*C, TransferFunction(1, 1, s))
    >>> F2.doit()
    TransferFunction(10*s**6 + 165*s**5 + 825*s**4 + 2005*s**3 + 2735*s**2 + 1880*s + 300, 11*s**6 + 189*s**5 + 1015*s**4 + 2617*s**3 + 3984*s**2 + 3260*s + 1200, s)

    To negate a ``Feedback`` object, the ``-`` operator can be prepended:

    >>> -F1
    Feedback(TransferFunction(-3*s**2 - 7*s + 3, s**2 - 4*s + 2, s), TransferFunction(5*s - 10, s + 7, s))
    >>> -F2
    Feedback(Series(TransferFunction(-1, 1, s), Series(TransferFunction(2*s**2 + 5*s + 1, s**2 + 2*s + 3, s), TransferFunction(5*s + 10, s + 10, s))), TransferFunction(1, 1, s))

    See Also
    ========

    TransferFunction, Series, Parallel

    """
    def __new__(cls, num, den):
        if not (isinstance(num, (TransferFunction, Series))
            and isinstance(den, (TransferFunction, Series))):
            raise TypeError("Unsupported type for numerator or denominator of Feedback.")

        if num == den:
            raise ValueError("The numerator cannot be equal to the denominator.")
        if not num.var == den.var:
            raise ValueError("Both numerator and denominator should be using the"
                " same complex variable.")
        obj = super().__new__(cls, num, den)
        obj._num = num
        obj._den = den
        obj._var = num.var

        return obj

    @property
    def num(self):
        """
        Returns the primary plant of the negative feedback closed loop.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> controller = TransferFunction(5*s - 10, s + 7, s)
        >>> F1 = Feedback(plant, controller)
        >>> F1.num
        TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> G = TransferFunction(2*s**2 + 5*s + 1, p**2 + 2*p + 3, p)
        >>> C = TransferFunction(5*p + 10, p + 10, p)
        >>> P = TransferFunction(1 - s, p + 2, p)
        >>> F2 = Feedback(TransferFunction(1, 1, p), G*C*P)
        >>> F2.num
        TransferFunction(1, 1, p)

        """
        return self._num

    @property
    def den(self):
        """
        Returns the feedback plant (often a feedback controller) of the
        negative feedback closed loop.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> controller = TransferFunction(5*s - 10, s + 7, s)
        >>> F1 = Feedback(plant, controller)
        >>> F1.den
        TransferFunction(5*s - 10, s + 7, s)
        >>> G = TransferFunction(2*s**2 + 5*s + 1, p**2 + 2*p + 3, p)
        >>> C = TransferFunction(5*p + 10, p + 10, p)
        >>> P = TransferFunction(1 - s, p + 2, p)
        >>> F2 = Feedback(TransferFunction(1, 1, p), G*C*P)
        >>> F2.den
        Series(TransferFunction(2*s**2 + 5*s + 1, p**2 + 2*p + 3, p), TransferFunction(5*p + 10, p + 10, p), TransferFunction(1 - s, p + 2, p))

        """
        return self._den

    @property
    def var(self):
        """
        Returns the complex variable of the Laplace transform used by all
        the transfer functions involved in the negative feedback closed loop.

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
        return self._var

    def doit(self, **kwargs):
        """
        Returns the resultant closed loop transfer function obtained by the
        negative feedback interconnection.

        Examples
        ========

        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> controller = TransferFunction(5*s - 10, s + 7, s)
        >>> F1 = Feedback(plant, controller)
        >>> F1.doit()
        TransferFunction(3*s**5 + 16*s**4 - 60*s**3 - 149*s**2 + 176*s - 42, 16*s**5 - 56*s**4 - 111*s**3 + 504*s**2 - 398*s + 88, s)
        >>> G = TransferFunction(2*s**2 + 5*s + 1, s**2 + 2*s + 3, s)
        >>> F2 = Feedback(G, TransferFunction(1, 1, s))
        >>> F2.doit()
        TransferFunction(2*s**4 + 9*s**3 + 17*s**2 + 17*s + 3, 3*s**4 + 13*s**3 + 27*s**2 + 29*s + 12, s)

        """
        arg_list = list(self.num.args) if isinstance(self.num, Series) else [self.num]
        # F_n and F_d are resultant TFs of num and den of Feedback.
        F_n, tf = self.num.doit(), TransferFunction(1, 1, self.num.var)
        F_d = Parallel(tf, Series(self.den, *arg_list)).doit()

        return TransferFunction(F_n.num*F_d.den, F_n.den*F_d.num, F_n.var)

    def _eval_rewrite_as_TransferFunction(self, num, den, **kwargs):
        return self.doit()

    def __neg__(self):
        return Feedback(-self.num, self.den)
