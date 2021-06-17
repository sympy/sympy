from sympy.matrices import ImmutableMatrix
from sympy import Basic, Mul, Pow, degree, Symbol, expand, cancel, Expr, exp, roots, ShapeError
from sympy.core.evalf import EvalfMixin
from sympy.core.logic import fuzzy_and
from sympy.core.numbers import Integer, ComplexInfinity
from sympy.core.sympify import sympify, _sympify
from sympy.polys import Poly, rootof
from sympy.series import limit

__all__ = ['TransferFunction', 'Series', 'Parallel', 'Feedback', 'TransferFunctionMatrix']


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

    Generally, a dynamical system respresenting a physical system can be described in terms of Linear
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

            $\small{b_{m}s^{m}\mathcal{L}[y]+b_{m-1}s^{m-1}\mathcal{L}[y]+\dots+b_{1}s\mathcal{L}[y]+b_{0}\mathcal{L}[y]=
            a_{n}s^{n}\mathcal{L}[x]+a_{n-1}s^{n-1}\mathcal{L}[x]+\dots+a_{1}s\mathcal{L}[x]+a_{0}\mathcal{L}[x]}$

    Here, the superscript on $s$ is exponent. Note that the zero initial conditions assumption, mentioned above, is very important
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
        number or a polynomial. Also, when ``num`` or ``den`` has a time
        delay term.
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

    >>> tf3 = TransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p)
    >>> tf3
    TransferFunction(p**2 + 2*p - 3, p**2 + 4*p - 5, p)

    To negate a transfer function the ``-`` operator can be prepended:

    >>> tf4 = TransferFunction(-a + s, p**2 + s, p)
    >>> -tf4
    TransferFunction(a - s, p**2 + s, p)
    >>> tf5 = TransferFunction(s**4 - 2*s**3 + 5*s + 4, s + 4, s)
    >>> -tf5
    TransferFunction(-s**4 + 2*s**3 - 5*s - 4, s + 4, s)

    You can use a Float or an Integer (or other constants) as numerator and denominator:

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
    TransferFunction(a**3 + 3*a**2*s + 3*a*s**2 + s**3, -a**3 + 3*a**2*s - 3*a*s**2 + s**3, s)
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
    Parallel(TransferFunction(-p + s, s + 3, s), Series(TransferFunction(-1, 1, s), Series(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(1 - s, s**2 + 4, s))))
    >>> tf11 * tf10 * tf9
    Series(TransferFunction(4*s**2 + 2*s - 4, s - 1, s), TransferFunction(-p + s, s + 3, s), TransferFunction(s + 1, s**2 + s + 1, s))
    >>> tf9 * tf11 + tf10 * tf12
    Parallel(Series(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(4*s**2 + 2*s - 4, s - 1, s)), Series(TransferFunction(-p + s, s + 3, s), TransferFunction(1 - s, s**2 + 4, s)))
    >>> (tf9 + tf12) * (tf10 + tf11)
    Series(Parallel(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(1 - s, s**2 + 4, s)), Parallel(TransferFunction(-p + s, s + 3, s), TransferFunction(4*s**2 + 2*s - 4, s - 1, s)))

    These unevaluated ``Series`` or ``Parallel`` objects can convert into the
    resultant transfer function using ``.doit()`` method or by ``.rewrite(TransferFunction)``.

    >>> ((tf9 + tf10) * tf12).doit()
    TransferFunction(p*s**3 - p - s**4 - s**3 - 3*s**2 + 2*s + 3, s**5 + 4*s**4 + 8*s**3 + 19*s**2 + 16*s + 12, s)
    >>> (tf9 * tf10 - tf11 * tf12).rewrite(TransferFunction)
    TransferFunction(-p*s**4 - 3*p*s**2 + 4*p + 4*s**6 + 15*s**5 + 2*s**4 - 13*s**3 - 14*s**2 - 6*s + 12, s**6 + 3*s**5 + 4*s**4 + 11*s**3 - 3*s**2 - 4*s - 12, s)

    See Also
    ========

    TransferFunctionMatrix, Feedback, Series, Parallel

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

        if (((isinstance(num, Expr) and not isinstance(num, ImmutableMatrix) and num.has(Symbol) and not num.has(exp)) or num.is_number) and
            ((isinstance(den, Expr) and not isinstance(num, ImmutableMatrix) and den.has(Symbol) and not den.has(exp)) or den.is_number)):
            obj = super(TransferFunction, cls).__new__(cls, num, den, var)
            obj._num = num
            obj._den = den
            obj._var = var
            obj._num_inputs, obj._num_outputs = 1, 1 # TODO: Remove _num_inputs and _num_outputs for TF and everywhere else except TFM
            return obj

        else:
            raise TypeError("Unsupported type for numerator or denominator of TransferFunction.")

    @classmethod
    def from_rational_expression(cls, expr, var=None):
        r"""
        Allows an alternative way for the users to instantiate a ``TransferFunction`` object.
        ``TransferFunction`` objects can be created simply by passing the rational expression.
        SymPy will smartly create a ``TransferFunction`` object with the properties of the
        expression.

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

            Wen ``expr`` has more than one variables and optional parameter
            ``var`` is not passed.
        ZeroDivisionError
            When denominator of ``expr`` is zero or it has ``ComplexInfinity``
            in its numerator.

        Examples
        ========

        >>> expr1 = (s + 5)/(3*s**2 + 2*s + 1)
        >>> tf1 = TransferFunction.from_rational_expression(expr1)
        >>> tf1
        TransferFunction(s + 5, 3*s**2 + 2*s + 1, s)
        >>> expr2 = (a*p**3 - a*p**2 + s*p)/(p + a**2)  # Expr with more than one variables
        >>> tf2 = TransferFunction(expr2, p)
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
        TransferFunction(a + a*s, s**2 + s + 1, s)

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
        if _den == 0 or _num.has(ComplexInfinity):
            raise ZeroDivisionError("TransferFunction can't have a zero denominator.")
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

    # TODO: Remove num_inputs, num_outputs and shape from TF

    # @property
    # def num_inputs(self):
    #     return self._num_inputs

    # @property
    # def num_outputs(self):
    #     return self._num_outputs

    # @property
    # def shape(self):
    #     return self._num_outputs, self._num_inputs

    # TODO: subs should return expr instead of TF

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

    # TODO: Implement expand(tf) by default and remove it from here.

    # def expand(self):
    #     """
    #     Returns the transfer function with numerator and denominator
    #     in expanded form.

    #     Examples
    #     ========

    #     >>> from sympy.abc import s, p, a, b
    #     >>> from sympy.physics.control.lti import TransferFunction
    #     >>> G1 = TransferFunction((a - s)**2, (s**2 + a)**2, s)
    #     >>> G1.expand()
    #     TransferFunction(a**2 - 2*a*s + s**2, a**2 + 2*a*s**2 + s**4, s)
    #     >>> G2 = TransferFunction((p + 3*b)*(p - b), (p - b)*(p + 2*b), p)
    #     >>> G2.expand()
    #     TransferFunction(-3*b**2 + 2*b*p + p**2, -2*b**2 + b*p + p**2, p)

    #     """
    #     return TransferFunction(expand(self.num), expand(self.den), self.var)

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
        elif other == -1:
            return -self  # Matrices compute -A internally by -1*element(A) for all element(A)
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

    # TODO: _to_expr() will convert TF to expr type. As Matrix addition/multiplication is only allowed for expr type.

    def _to_expr(self):
        """
        To convert TransferFunction type to SymPy Expr
        """
        return Mul(self.num, Pow(self.den, -1, evaluate=False), evaluate=False)


class Series(Basic):
    r"""
    A class for representing series configuration of SISO and MIMO transfer functions.

    Parameters
    ==========

    args : TransferFunction, TransferFunctionMatrix, Series, Parallel
        Systems in series configuration.
    evaluate : Boolean, Keyword
        When passed ``True``, returns the equivalent
        Series(*args).doit(). Set to ``False`` by default.

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
        type of systems passed.

    Examples
    ========

    SISO-System Examples -

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

    MIMO-System Examples-

    Notes
    =====

    All the transfer functions should use the same complex variable ``var`` of the Laplace transform.

    ``Series(A, B)`` is not equivalent to ``A*B``. It is always in the reverse order, that is ``B*A``.
    In, SISO systems, it hardly matters, but in MIMO systems, there can be issues due to non-commutative
    nature of Matrix Multiplication.

    See Also
    ========

    Parallel, Feedback

    """

    # TODO: Series implementation for TFM

    def __new__(cls, *args, evaluate=False):
        # if len(args) == 0:
        #     raise ValueError("Needs at least 1 argument.")
        # if not all(isinstance(arg, (TransferFunction, TransferFunctionMatrix, Parallel, Series))
        #     for arg in args):
        #     raise TypeError("Unsupported type of argument(s) for Series.")

        # obj = super(Series, cls).__new__(cls, *args)
        # obj._is_not_matrix = all(isinstance(arg.doit(), TransferFunction) for arg in args)
        # if not obj._is_not_matrix:
        #     obj._num_outputs, obj._num_inputs = args[0].num_outputs, args[-1].num_inputs
        #     for x in range(len(args) - 1):
        #         # input-output sizes should be consistent.
        #         if args[x].num_inputs != args[x + 1].num_outputs:
        #             raise ValueError("Argument {0} of Series has {1} input(s),"
        #                 " but argument {2} has {3} output(s)."
        #                 .format(x + 1, args[x].num_inputs, x + 2, args[x + 1].num_outputs))
        # else:
        #     obj._num_outputs, obj._num_inputs = 1, 1

        # tf = "transfer functions" if obj._is_not_matrix else "TransferFunctionMatrix objects"
        # obj._var = args[0].var
        # if not all(arg.var == obj._var for arg in args):
        #     raise ValueError("All {0} should use the same complex"
        #         " variable of the Laplace transform.".format(tf))
        # return obj.doit() if evaluate else obj

        if len(args) == 0:
            raise ValueError("Needs at least 1 argument.")
        if not all(isinstance(arg, (TransferFunction, TransferFunctionMatrix, Parallel, Series)) for arg in args):
            raise TypeError("Unsupported type of argument(s) for Series.")

        if all(isinstance(arg, TransferFunction) for arg in args):
            obj = super().__new__(cls, *args)
            obj._var = None
            obj._type_func = "TransferFunction"
            for arg in args:
                if obj._var is None:
                    obj._var = arg.var
                elif obj._var != arg.var:
                    raise ValueError("All transfer functions should use the same complex"
                        " variable of the Laplace transform.")
            if evaluate:
                return obj.doit()
            return obj
        elif all(isinstance(arg, TransferFunctionMatrix) for arg in args):
            obj = super().__new__(cls, *args)
            obj._var = None
            obj._type_func = "TransferFunctionMatrix"

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

    @property
    def num_inputs(self):
        # If the Series is for TFMs, then return num_inputs of the first TFM arg as num_input for the series system
        try:
            return self.args[0].num_inputs
        # If no such attribute is found, then its TF instead of TFM and num inputs for tf should be None
        except AttributeError:
            return None

    @property
    def num_outputs(self):
        try:
            return self.args[len(self.args) - 1].num_outputs
        except AttributeError:
            return None

    @property
    def shape(self):
        try:
            return self.args[len(self.args) - 1].num_outputs, self.args[0].num_inputs
        finally:
            return None

    # TODO: Implement Series().doit() for tfm. Matrix multiplication of ImmutableMatrices can be used

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
        TransferFunction((s**3 - 2)*(a*p**2 + b*s), (-p + s)*(s**4 + 5*s + 6), s)
        >>> Series(-tf1, -tf2).doit()
        TransferFunction((2 - s**3)*(-a*p**2 - b*s), (-p + s)*(s**4 + 5*s + 6), s)

        """
        res = None
        for arg in self.args:
            arg = arg.doit()
            if res is None:
                res = arg
            else:
                if self._is_not_matrix:
                    num_ = arg.num * res.num
                    den_ = arg.den * res.den
                    res = TransferFunction(num_, den_, self.var)
                else:
                    # Now we multiply two transfer function matrices...
                    # If we have two TFMs with shape (a, b) and (b, d), respectively,
                    # then the resultant TFM will be of shape (a, d).
                    if res.num_outputs == 1 and arg.num_inputs == 1:
                        a = [None]
                        for i in range(res.num_inputs):
                            if a[0] is None:
                                a[0] = res.args[0][0][i] * arg.args[0][i]
                            else:
                                a[0] += res.args[0][0][i] * arg.args[0][i]

                    elif res.num_outputs == 1 or arg.num_inputs == 1:
                        if res.num_outputs == 1:
                            a = [[None] * arg.num_inputs]

                            for j in range(arg.num_inputs):
                                for k in range(arg.num_outputs):
                                    if a[0][j] is None:
                                        a[0][j] = res.args[0][0][j] * arg.args[0][k][j]
                                    else:
                                        a[0][j] += res.args[0][0][j] * arg.args[0][k][j]

                        else:
                            a = [None] * res.num_outputs

                            for i in range(res.num_outputs):
                                for k in range(arg.num_outputs):
                                    if a[i] is None:
                                        if res.num_inputs == 1:
                                            a[i] = res.args[0][i] * arg.args[0][k]
                                        else:
                                            a[i] = res.args[0][i][k] * arg.args[0][k]
                                    else:
                                        a[i] += res.args[0][i][k] * arg.args[0][k]

                    else:
                        a = [[None] * arg.num_inputs for _ in range(res.num_outputs)]
                        for i in range(res.num_outputs):
                            for j in range(arg.num_inputs):
                                for k in range(arg.num_outputs):
                                    if a[i][j] is None:     # First operation.
                                        if res.num_inputs == 1:
                                            a[i][j] = res.args[0][i] * arg.args[0][k][j]
                                        else:
                                            a[i][j] = res.args[0][i][k] * arg.args[0][k][j]
                                    else:
                                        a[i][j] += res.args[0][i][k] * arg.args[0][k][j]

                    res = TransferFunctionMatrix(a, (res.num_outputs, arg.num_inputs), arg.var)

        return res

    def _eval_rewrite_as_TransferFunction(self, *args, **kwargs):
        """ Series(tfm1, tfm2, Parallel(tfm3, tfm4, ...), ...) not allowed. """
        if not self._is_not_matrix:
            raise ValueError("Only transfer functions or a collection of transfer functions"
                " is allowed in the arguments.")

        return self.doit()

    def _eval_rewrite_as_TransferFunctionMatrix(self, *args, **kwargs):
        """ Series(tf1, tf2, Parallel(tf3, tf4, ...), ...) not allowed. """
        if self._is_not_matrix:
            raise ValueError("Only transfer function matrices or a collection of transfer function"
                " matrices is allowed in the arguments.")

        return self.doit()

    def __add__(self, other):
        if isinstance(other, (TransferFunction, Series, TransferFunctionMatrix)):
            if isinstance(other, Series):
                if self._is_not_matrix != other._is_not_matrix:
                    raise ValueError("Both Series objects should either handle SISO or MIMO"
                        " transfer function.")

            if isinstance(other, TransferFunctionMatrix):
                if self._is_not_matrix:
                    raise ValueError("Series object should handle MIMO transfer function to"
                        " perform this addition.")
            if not self.shape == other.shape:
                raise ShapeError("Shapes of operands are not compatible for addition.")
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")

            return Parallel(self, other)
        elif isinstance(other, Parallel):
            if self._is_not_matrix != other._is_not_matrix:
                raise ValueError("Both Series and Parallel objects should either handle SISO or MIMO"
                    " transfer function.")

            if not self.shape == other.shape:
                raise ShapeError("Shapes of operands are not compatible for addition.")
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            arg_list = list(other.args)

            return Parallel(self, *arg_list)
        else:
            raise ValueError("This expression is invalid.")

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if isinstance(other, (TransferFunction, Parallel, TransferFunctionMatrix)):
            if isinstance(other, Parallel):
                if self._is_not_matrix != other._is_not_matrix:
                    raise ValueError("Both Series and Parallel objects should either handle SISO or MIMO"
                        " transfer function.")
            if isinstance(other, TransferFunctionMatrix):
                if self._is_not_matrix:
                    raise ValueError("Series object should only handle MIMO transfer function to"
                        " perform this multiplication.")

            if not self.num_inputs == other.num_outputs:
                raise ValueError("C = A * B: A has {0} input(s), but B has {1} output(s)."
                    .format(self.num_inputs, other.num_outputs))
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            arg_list = list(self.args)

            return Series(*arg_list, other)
        elif isinstance(other, Series):
            if self._is_not_matrix != other._is_not_matrix:
                raise ValueError("Both Series objects should either handle SISO or MIMO"
                    " transfer function.")

            if not self.num_inputs == other.num_outputs:
                raise ValueError("C = A * B: A has {0} input(s), but B has {1} output(s)."
                    .format(self.num_inputs, other.num_outputs))
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")
            self_arg_list = list(self.args)
            other_arg_list = list(other.args)

            return Series(*self_arg_list, *other_arg_list)
        else:
            raise ValueError("This expression is invalid.")

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
            raise ValueError("This expression is invalid.")

    def __neg__(self):
        if self._is_not_matrix:
            return Series(TransferFunction(-1, 1, self.var), self)
        else:
            neg_tfm = [[TransferFunction(-1, 1, self.var)] * self.num_outputs]
            return Series(TransferFunctionMatrix(neg_tfm), self)

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

    @property
    def is_SISO(self):
        if self._is_not_matrix:
            return True
        else:
            return False


class Parallel(Basic):
    r"""
    A class for representing parallel configuration of SISO and MIMO transfer functions.

    Parameters
    ==========

    args : TransferFunction, TransferFunctionMatrix, Series, Parallel
        Systems in parallel arrangement
    evaluate : Boolean, Keyword
        When passed ``True``, returns the equivalent
        Parallel(*args).doit(). Set to ``False`` by default.

    Raises
    ======

    ValueError
        When no argument is passed.

        ``var`` attribute is not same for every system.

        All MIMO systems passed don't have same shape.
    TypeError
        Any of the passed ``*args`` has unsupported type

        A combination of SISO and MIMO systems is
        passed. There should be homogeneity in the
        type of systems passed.

    Examples
    ========

    SISO-System examples -

    >>> from sympy.abc import s, p, a, b
    >>> from sympy.printing import pprint
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
    TransferFunction(-p**2*(-p + s)*(s**4 + 5*s + 6) + (p + s)*((-p + s)*(s**3 - 2) + (a*p**2 + b*s)*(s**4 + 5*s + 6)), (-p + s)*(p + s)*(s**4 + 5*s + 6), s)
    >>> Parallel(tf2, Series(tf1, -tf3)).doit()
    TransferFunction(-p**2*(a*p**2 + b*s)*(s**4 + 5*s + 6) + (-p + s)*(p + s)*(s**3 - 2), (-p + s)*(p + s)*(s**4 + 5*s + 6), s)

    MIMO-System example -

    >>> expr_1 = 1/s
    >>> expr_2 = s/(s**2-1)
    >>> expr_3 = (2 + s)/(s**2 - 1)
    >>> expr_4 = 5
    >>> tfm_a = TransferFunctionMatrix([[expr_1, expr_2], [expr_3, expr_4]])
    >>> tfm_b = TransferFunctionMatrix([[expr_2, expr_1], [expr_4, expr_3]])
    >>> tfm_c = TransferFunctionMatrix([[expr_3, expr_4], [expr_1, expr_2]])
    >>> Parallel(tfm_a, tfm_b, tfm_c)
    Parallel(Matrix([
    [           TransferFunction(1, s, s), TransferFunction(s, s**2 - 1, s)],
    [TransferFunction(s + 2, s**2 - 1, s),        TransferFunction(5, 1, s)]]), Matrix([
    [TransferFunction(s, s**2 - 1, s),            TransferFunction(1, s, s)],
    [       TransferFunction(5, 1, s), TransferFunction(s + 2, s**2 - 1, s)]]), Matrix([
    [TransferFunction(s + 2, s**2 - 1, s),        TransferFunction(5, 1, s)],
    [           TransferFunction(1, s, s), TransferFunction(s, s**2 - 1, s)]]))
    >>> pprint(_, use_unicode=False)  #  For Better Visualization
    [  1       s   ]   [  s       1   ]   [s + 2     5   ]
    [  -     ------]   [------    -   ]   [------    -   ]
    [  s      2    ]   [ 2        s   ]   [ 2        1   ]
    [        s  - 1]   [s  - 1        ]   [s  - 1        ]
    [              ] + [              ] + [              ]
    [s + 2     5   ]   [  5     s + 2 ]   [  1       s   ]
    [------    -   ]   [  -     ------]   [  -     ------]
    [ 2        1   ]   [  1      2    ]   [  s      2    ]
    [s  - 1        ]   [        s  - 1]   [        s  - 1]
    >>> Parallel(tfm_a, tfm_b, tfm_c).doit()
    Matrix([[TransferFunction(3*s**2 + 2*s - 1, s**3 - s, s), TransferFunction(5*s**3 + 2*s**2 - 5*s - 1, s**3 - s, s)], [TransferFunction(5*s**3 + 2*s**2 - 3*s - 1, s**3 - s, s), TransferFunction(5*s**2 + 2*s - 3, s**2 - 1, s)]])
    >>> pprint(_, use_unicode=False)
    [      2                   3      2          ]
    [   3*s  + 2*s - 1      5*s  + 2*s  - 5*s - 1]
    [   --------------      ---------------------]
    [        3                       3           ]
    [       s  - s                  s  - s       ]
    [                                            ]
    [   3      2                  2              ]
    [5*s  + 2*s  - 3*s - 1     5*s  + 2*s - 3    ]
    [---------------------     --------------    ]
    [         3                     2            ]
    [        s  - s                s  - 1        ]

    Notes
    =====

    All the transfer functions should use the same complex variable
    ``var`` of the Laplace transform.

    See Also
    ========

    Series, TransferFunction, Feedback

    """
    def __new__(cls, *args, evaluate=False):
        if len(args) == 0:
            raise ValueError("Needs at least 1 argument.")
        if not all(isinstance(arg, (TransferFunction, TransferFunctionMatrix, Series, Parallel))
            for arg in args):
            raise TypeError("Unsupported type of argument(s) for Parallel.")
        if not all(arg.var == args[0].var for arg in args):
            raise ValueError("`var` attribute of each arg should be same")

        def _instantiate_SISO_Parallel(*args):
            obj = super(Parallel, cls).__new__(cls, *args)
            obj._var = args[0].var
            obj._is_SISO = True
            obj._shape = None
            obj._num_inputs, obj._num_outputs = None, None
            return obj.doit() if evaluate else obj

        def _instantiate_MIMO_Parallel(*args):
            if not all(args[0].shape == arg.shape for arg in args):
                raise ValueError("All MIMO arguments should have same shape.")
            obj = super(Parallel, cls).__new__(cls, *args)
            obj._var = args[0].var
            obj._is_SISO = False
            obj._shape = args[0].shape
            obj._num_inputs, obj._num_outputs = args[0].shape
            return obj.doit() if evaluate else obj

        _num_SISO, _num_MIMO = 0, 0
        for arg in args:
            if isinstance(arg, TransferFunction):
                _num_SISO += 1
            elif isinstance(arg, TransferFunctionMatrix):
                _num_MIMO += 1
            elif arg.is_SISO:
                _num_SISO += 1
            else:
                _num_MIMO += 1

        if _num_MIMO == 0 or _num_SISO == 0:
            if _num_MIMO == 0:
                return _instantiate_SISO_Parallel(*args)
            if _num_SISO == 0:
                return _instantiate_MIMO_Parallel(*args)
        else:
            raise TypeError("Not a valid `Parallel` object. Either all the arguments should be"
                " SISO or MIMO, not a combination of both.")

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

    @property
    def is_SISO(self):
        return True if self._is_SISO else False

    @property
    def num_inputs(self):
        return self._num_inputs

    @property
    def num_outputs(self):
        return self._num_outputs

    @property
    def shape(self):
        return self._shape

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
        TransferFunction((-p + s)*(s**3 - 2) + (a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s)
        >>> Parallel(-tf1, -tf2).doit()
        TransferFunction((2 - s**3)*(-p + s) + (-a*p**2 - b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s)

        """
        _args = list(self.args)
        def _SISO_doit():
            res = 0
            for arg_index in range(len(_args)):
                if not isinstance(_args[arg_index], TransferFunction):
                    _args[arg_index] = _args[arg_index].doit()

                res += _args[arg_index]._to_expr()
            return TransferFunction.from_rational_expression(res, self.var)

        def _MIMO_doit():
            res = _args[0].doit()._to_Immutable_Matrix(expr=True)
            for arg_index in range(len(_args) - 1):
                if not isinstance(_args[arg_index + 1], TransferFunctionMatrix):
                    _args[arg_index + 1] = _args[arg_index + 1].doit()

                res += _args[arg_index + 1]._to_Immutable_Matrix(expr=True)
            return TransferFunctionMatrix(*res.args, var=self.var)

        if self.is_SISO:
            return _SISO_doit()
        else:
            return _MIMO_doit()

    def _eval_rewrite_as_TransferFunction(self, *args, **kwargs):
        """ Parallel(tfm1, tfm2, Series(tfm3, tfm4, ...), ...) not allowed. """
        if not self.is_SISO:
            raise ValueError("Only transfer functions or a collection of transfer functions"
                " is allowed in the arguments.")

        return self.doit()

    def _eval_rewrite_as_TransferFunctionMatrix(self, *args, **kwargs):
        """ Parallel(tf1, tf2, Series(tf3, tf4, ...), ...) not allowed. """
        if self.is_SISO:
            raise ValueError("Only transfer function matrices or a collection of transfer function"
                " matrices is allowed in the arguments.")

        return self.doit()

    def __add__(self, other):
        if isinstance(other, (TransferFunction, Series, TransferFunctionMatrix)):
            arg_list = list(self.args)
            arg_list.append(other)
            return Parallel(*arg_list)
        elif isinstance(other, Parallel):
            self_arg_list = list(self.args)
            other_arg_list = list(other.args)
            for elem in other_arg_list:
                self_arg_list.append(elem)
            return Parallel(*self_arg_list)
        else:
            raise ValueError("This expression is invalid.")

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, (TransferFunction, Parallel, TransferFunctionMatrix)):
            if isinstance(other, Parallel):
                if self._is_not_matrix != other._is_not_matrix:
                    raise ValueError("Both Parallel objects should either handle SISO or MIMO"
                        " transfer function.")
            if isinstance(other, TransferFunctionMatrix):
                if self._is_not_matrix:
                    raise ValueError("The Parallel object should only handle MIMO transfer function to"
                        " perform this multiplication.")
            if not self.num_inputs == other.num_outputs:
                raise ValueError("C = A * B: A has {0} input(s), but B has {1} output(s)."
                    .format(self.num_inputs, other.num_outputs))
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")

            return Series(self, other)
        elif isinstance(other, Series):
            if self._is_not_matrix != other._is_not_matrix:
                raise ValueError("Both Series and Parallel objects should either handle SISO or MIMO"
                    " transfer function.")
            if not self.num_inputs == other.num_outputs:
                raise ValueError("C = A * B: A has {0} input(s), but B has {1} output(s)."
                    .format(self.num_inputs, other.num_outputs))
            if not self.var == other.var:
                raise ValueError("All the transfer functions should use the same complex variable "
                    "of the Laplace transform.")

            arg_list = list(other.args)
            return Series(self, *arg_list)
        else:
            raise ValueError("This expression is invalid.")

    __rmul__ = __mul__

    def __neg__(self):
        if self._is_not_matrix:
            return Series(TransferFunction(-1, 1, self.var), self)
        else:
            neg_args = [-arg for arg in self.args]
            return Parallel(*neg_args)

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
    TransferFunction((s + 7)*(s**2 - 4*s + 2)*(3*s**2 + 7*s - 3), ((s + 7)*(s**2 - 4*s + 2) + (5*s - 10)*(3*s**2 + 7*s - 3))*(s**2 - 4*s + 2), s)
    >>> G = TransferFunction(2*s**2 + 5*s + 1, s**2 + 2*s + 3, s)
    >>> C = TransferFunction(5*s + 10, s + 10, s)
    >>> F2 = Feedback(G*C, TransferFunction(1, 1, s))
    >>> F2.doit()
    TransferFunction((s + 10)*(5*s + 10)*(s**2 + 2*s + 3)*(2*s**2 + 5*s + 1), (s + 10)*((s + 10)*(s**2 + 2*s + 3) + (5*s + 10)*(2*s**2 + 5*s + 1))*(s**2 + 2*s + 3), s)

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
        TransferFunction((s + 7)*(s**2 - 4*s + 2)*(3*s**2 + 7*s - 3), ((s + 7)*(s**2 - 4*s + 2) + (5*s - 10)*(3*s**2 + 7*s - 3))*(s**2 - 4*s + 2), s)
        >>> G = TransferFunction(2*s**2 + 5*s + 1, s**2 + 2*s + 3, s)
        >>> F2 = Feedback(G, TransferFunction(1, 1, s))
        >>> F2.doit()
        TransferFunction((s**2 + 2*s + 3)*(2*s**2 + 5*s + 1), (s**2 + 2*s + 3)*(3*s**2 + 7*s + 4), s)

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


class TransferFunctionMatrix(ImmutableMatrix, Basic):
    r"""
    A class for representing the MIMO (multiple-input and multiple-output)
    generalization of the SISO (single-input and single-output) transfer function.

    It is a matrix of transfer functions (``TransferFunction`` objects).
    The arguments are ``arg`` and ``var``, where ``arg``
    is the compulsory argument. ``arg`` is expected to be of the type list or list of lists
    which holds the transfer functions. ``var``
    is a complex variable (which has to be a ``Symbol``) of the Laplace transform
    used by all the transfer functions in the matrix.

    Parameters
    ==========

    arg: List, Nested List
        Users are expected to input a nested list of ``TransferFunction``
        (or ``Expr``) objects like they input Numbers in a normal SymPy matrix.
        In case a user inputs a List, it would be considered as a column matrix
        and not row. So, its better to pass nested list to avoid such dilemma.
    var: Symbol, keyword, optional
        ``var`` is an optional keyword argument. If ``var`` is not passed explicitly
        by the user then ``var`` of the first element is taken by default.

    Examples
    ========

    ``pprint()`` can be used for better visualization of ``TransferFunctionMatrix`` objects.

    >>> from sympy.abc import s, p, a
    >>> from sympy.printing import pprint
    >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
    >>> tf_1 = TransferFunction(s + a, s**2 + s + 1, s)
    >>> tf_2 = TransferFunction(p**4 - 3*p + 2, s + p, s)
    >>> tf_3 = TransferFunction(3, s + 2, s)
    >>> tf_4 = TransferFunction(-a + p, 9*s - 9, s)
    >>> tfm_1 = TransferFunctionMatrix([tf_1, tf_2, tf_3])
    >>> tfm_1
    Matrix([[TransferFunction(a + s, s**2 + s + 1, s)], [TransferFunction(p**4 - 3*p + 2, p + s, s)], [TransferFunction(3, s + 2, s)]])
    >>> tfm_1.var
    s
    >>> tfm_1.num_inputs
    1
    >>> tfm_1.num_outputs
    3
    >>> tfm_1.shape
    (3, 1)
    >>> tfm_1.args  # Structurally similar to ImmutableDenseMatrix().args
    (3, 1, (TransferFunction(a + s, s**2 + s + 1, s), TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(3, s + 2, s)))
    >>> tfm_2 = TransferFunctionMatrix([[tf_1, -tf_3], [tf_2, -tf_1], [tf_3, -tf_2]])
    >>> tfm_2
    Matrix([[TransferFunction(a + s, s**2 + s + 1, s), TransferFunction(-3, s + 2, s)], [TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(-a - s, s**2 + s + 1, s)], [TransferFunction(3, s + 2, s), TransferFunction(-p**4 + 3*p - 2, p + s, s)]])
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
    [   s + 2          p + s     ]

    TransferFunctionMatrix can be transposed, if user wants to switch the input and output transfer functions

    >>> tfm_2.transpose()
    Matrix([[TransferFunction(a + s, s**2 + s + 1, s), TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(3, s + 2, s)], [TransferFunction(-3, s + 2, s), TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(-p**4 + 3*p - 2, p + s, s)]])
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
    [             s  + s + 1                 ]

    TransferFunctionMatrix objects can also be instantiated by passing SymPy expressions

    >>> expr_1 = 5/s
    >>> expr_2 = 5*s/(2 + s**2)
    >>> expr_3 = 5/((s)*(2 + s**2))
    >>> expr_4 = 5
    >>> tfm_3 = TransferFunctionMatrix([[expr_1, expr_2], [expr_3, expr_4]])  # SymPy will assume the var of the first element as the var for the system
    >>> tfm_3
    Matrix([[TransferFunction(5, s, s), TransferFunction(5*s, s**2 + 2, s)], [TransferFunction(5, s**3 + 2*s, s), TransferFunction(5, 1, s)]])
    >>> pprint(tfm_3, use_unicode=False)
    [   5       5*s  ]
    [   -      ------]
    [   s       2    ]
    [          s  + 2]
    [                ]
    [   5        5   ]
    [--------    -   ]
    [ 3          1   ]
    [s  + 2*s        ]
    >>> tfm_3.var
    s
    >>> tfm_3.shape
    (2, 2)
    >>> tfm_3.num_outputs
    2
    >>> tfm_3.num_inputs
    2
    >>> tfm_3.args
    (2, 2, (TransferFunction(5, s, s), TransferFunction(5*s, s**2 + 2, s), TransferFunction(5, s**3 + 2*s, s), TransferFunction(5, 1, s)))

    To access the ``TransferFunction`` at any index in the ``TransferFunctionMatrix``, use the index notation.

    >>> tfm_3[1, 0]  # gives the TransferFunction present at 2nd Row and 1st Col. Similar to that in Matrix classes
    TransferFunction(5, s**3 + 2*s, s)
    >>> tfm_3[0, 0]  # gives the TransferFunction present at 1st Row and 1st Col.
    TransferFunction(5, s, s)
    >>> tfm_3[:, 0]  # gives the first column
    Matrix([[TransferFunction(5, s, s)], [TransferFunction(5, s**3 + 2*s, s)]])
    >>> pprint(_, use_unicode=False)
    [   5    ]
    [   -    ]
    [   s    ]
    [        ]
    [   5    ]
    [--------]
    [ 3      ]
    [s  + 2*s]
    >>> tfm_3[0, :]  # gives the first row
    Matrix([[TransferFunction(5, s, s), TransferFunction(5*s, s**2 + 2, s)]])
    >>> pprint(_, use_unicode=False)
    [5   5*s  ]
    [-  ------]
    [s   2    ]
    [   s  + 2]

    To negate a transfer function matrix, ``-`` operator can be prepended:

    >>> tfm_4 = TransferFunctionMatrix([tf_2, -tf_1, tf_3])  # Remember-Simple lists are column matrix
    >>> -tfm_4
    Matrix([[TransferFunction(-p**4 + 3*p - 2, p + s, s)], [TransferFunction(a + s, s**2 + s + 1, s)], [TransferFunction(-3, s + 2, s)]])
    >>> tfm_5 = TransferFunctionMatrix([[tf_1, tf_2], [tf_3, -tf_1]])
    >>> -tfm_5
    Matrix([[TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(-p**4 + 3*p - 2, p + s, s)], [TransferFunction(-3, s + 2, s), TransferFunction(a + s, s**2 + s + 1, s)]])

    ``subs()`` returns an ImmutableDenseMatrix object with the value substituted in the expression. This will not
    mutate your original ``TransferFunctionMatrix``.

    >>> from sympy import I
    >>> tfm_3.subs({s: I})  #  substituting imaginary i (iota) to s
    Matrix([[-5*I, 5*I], [-5*I, 5/1]])
    >>> type(_)
    <class 'sympy.matrices.immutable.ImmutableDenseMatrix'>
    >>> type(tf_3)
    <class 'sympy.physics.control.lti.TransferFunctionMatrix'>

    Addition, subtraction, and multiplication of transfer function matrices can form
    unevaluated ``Series`` or ``Parallel`` objects.

    - For addition and subtraction:
      All the transfer function matrices must have the same shape.

    - For multiplication (C = A * B):
      The number of inputs of the first transfer function matrix (A) must be equal to the
      number of outputs of the second transfer function matrix (B).

    >>> tfm_6 = TransferFunctionMatrix([[tf_4, -tf_1, tf_3], [-tf_2, -tf_4, -tf_3]])
    >>> tfm_7 = TransferFunctionMatrix([tf_3, tf_2, -tf_1])
    >>> tfm_8 = TransferFunctionMatrix([-tf_3])
    >>> tfm_9 = TransferFunctionMatrix([tf_1, tf_2, tf_4])
    >>> tfm_10 = TransferFunctionMatrix([tf_4, -tf_1])
    >>> tfm_7 + tfm_9
    Parallel(TransferFunctionMatrix([TransferFunction(3, s + 2, s), TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(-a - s, s**2 + s + 1, s)]), TransferFunctionMatrix((TransferFunction(a + s, s**2 + s + 1, s), TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(-a + p, 9*s - 9, s))))
    >>> -tfm_10 - tfm_8
    Parallel(TransferFunctionMatrix([TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(-p**4 + 3*p - 2, p + s, s), TransferFunction(a - p, 9*s - 9, s)]), TransferFunctionMatrix([TransferFunction(-3, s + 2, s), TransferFunction(-p**4 + 3*p - 2, p + s, s), TransferFunction(a + s, s**2 + s + 1, s)]))
    >>> tfm_7 * tfm_8
    Series(TransferFunctionMatrix([[TransferFunction(-a + p, 9*s - 9, s), TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(3, s + 2, s)], [TransferFunction(-p**4 + 3*p - 2, p + s, s), TransferFunction(a - p, 9*s - 9, s), TransferFunction(-3, s + 2, s)]]), TransferFunctionMatrix([TransferFunction(3, s + 2, s), TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(-a - s, s**2 + s + 1, s)]))
    >>> tfm_7 * tfm_8 * tfm_9
    Series(TransferFunctionMatrix([[TransferFunction(-a + p, 9*s - 9, s), TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(3, s + 2, s)], [TransferFunction(-p**4 + 3*p - 2, p + s, s), TransferFunction(a - p, 9*s - 9, s), TransferFunction(-3, s + 2, s)]]), TransferFunctionMatrix([TransferFunction(3, s + 2, s), TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(-a - s, s**2 + s + 1, s)]), TransferFunctionMatrix([TransferFunction(-3, s + 2, s)]))
    >>> tfm_10 + tfm_8*tfm_9
    Parallel(TransferFunctionMatrix((TransferFunction(a + s, s**2 + s + 1, s), TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(-a + p, 9*s - 9, s))), Series(TransferFunctionMatrix([TransferFunction(3, s + 2, s), TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(-a - s, s**2 + s + 1, s)]), TransferFunctionMatrix([TransferFunction(-3, s + 2, s)])))

    These unevaluated ``Series`` or ``Parallel`` objects can convert into the
    resultant transfer function matrix using ``.doit()`` method or by
    ``.rewrite(TransferFunctionMatrix)``.

    >>> (-tfm_8 + tfm_10 + tfm_8*tfm_9).doit()
    TransferFunctionMatrix([Parallel(TransferFunction(-3, s + 2, s), TransferFunction(a + s, s**2 + s + 1, s), Series(TransferFunction(3, s + 2, s), TransferFunction(-3, s + 2, s))), Parallel(TransferFunction(-p**4 + 3*p - 2, p + s, s), TransferFunction(p**4 - 3*p + 2, p + s, s), Series(TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(-3, s + 2, s))), Parallel(TransferFunction(a + s, s**2 + s + 1, s), TransferFunction(-a + p, 9*s - 9, s), Series(TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(-3, s + 2, s)))])
    >>> (-tfm_7 * -tfm_8 * -tfm_9).rewrite(TransferFunctionMatrix)
    TransferFunctionMatrix([Series(Parallel(Series(TransferFunction(a - p, 9*s - 9, s), TransferFunction(-3, s + 2, s)), Series(TransferFunction(a + s, s**2 + s + 1, s), TransferFunction(-p**4 + 3*p - 2, p + s, s)), Series(TransferFunction(-3, s + 2, s), TransferFunction(a + s, s**2 + s + 1, s))), TransferFunction(3, s + 2, s)), Series(Parallel(Series(TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(-3, s + 2, s)), Series(TransferFunction(-a + p, 9*s - 9, s), TransferFunction(-p**4 + 3*p - 2, p + s, s)), Series(TransferFunction(3, s + 2, s), TransferFunction(a + s, s**2 + s + 1, s))), TransferFunction(3, s + 2, s))])

    See Also
    ========

    TransferFunction, Series, Parallel, Feedback

    """
    def __new__(cls, *arg, var=None):
        # if not (isinstance(arg, (list, tuple)) and
        #     (all(isinstance(i, (list, tuple)) for i in arg) or
        #     all(isinstance(i, (TransferFunction, Parallel, Series, Expr, int, float)) for i in arg))):
        #     raise TypeError("Unsupported type for argument of TransferFunctionMatrix.")

        if len(arg) == 0:
            raise ValueError("Positional argument `arg` not passed.")

        # Converting flat tuple to nested list.
        if len(arg) > 1:
            def group(flat, size): return [flat[i:i+size] for i in range(0, len(flat), size)]
            tfm_cols = arg[1]
            flat_list = list(arg[2])
            arg = group(flat_list, tfm_cols)

        if isinstance(arg, tuple):
            arg = arg[0]

        # Checking if the list is NOT nested and then converting to nested if not already
        if not any(isinstance(i, list) for i in arg):
            arg = [[i] for i in arg]

        if not (isinstance(arg, list)):
            raise TypeError("Unsupported type for argument of TransferFunctionMatrix.")

        if not all(isinstance(arg[row][col], (TransferFunction, Series, Parallel, Expr, int, float))
                for col in range(len(arg[0])) for row in range(len(arg))):
                raise TypeError("All the lists in the first argument of TransferFunctionMatrix"
                    " only support transfer functions, and Series/Parallel objects in them.")
        if not all(len(l) == len(arg[0]) for l in arg):
                raise ValueError("Length of all the lists in the argument of"
                    " TransferFunctionMatrix should be equal.")

        if not var:
            if not isinstance(arg[0][0], (TransferFunction, Series, Parallel)):
                try:
                    var = TransferFunction.from_rational_expression(arg[0][0]).var
                except ValueError:
                    raise ValueError("Input the parameter `var` explicitly for this TransferFunction.")
            else:
                var = arg[0][0].var

        if not isinstance(var, Symbol):
            raise TypeError("`var` must be a Symbol, not {}.".format(type(var)))

        for row_num in range(len(arg)):
            for col_num in range(len(arg[row_num])):
                if not isinstance(arg[row_num][col_num], (TransferFunction, Series, Parallel)):
                    arg[row_num][col_num] = TransferFunction.from_rational_expression(arg[row_num][col_num], var=var)
                else:
                    if (arg[row_num][col_num]).var != var:
                        raise ValueError("Conflicting value(s) found for `var`. All TransferFunction instances in "
                            "TransferFunctionMatrix should use the same complex variable in Laplace domain.")
                    if isinstance(arg[row_num][col_num], (Series, Parallel)) and not arg[row_num][col_num].is_SISO:
                            raise TypeError("Only SISO type `Series`/`Parallel` are allowed in TransferFunctionMatrix")


        obj = super(TransferFunctionMatrix, cls).__new__(cls, arg)
        obj._var = arg[0][0].var
        obj._num_outputs, obj._num_inputs = obj.shape
        return obj

    @property
    def var(self):
        """
        Returns the complex variable used by all the transfer functions or
        ``Series``/``Parallel`` objects in a transfer function matrix.

        If explicitly provided as the third argument of ``TransferFunctionMatrix``, then returns
        the same, else, ``TransferFunctionMatrix`` automatically looks up the ``var`` used by
        transfer functions inside the list/tuple provided as the first argument and returns that.

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
        >>> tfm1 = TransferFunctionMatrix([G1, G2, G3], (3, 1), p)
        >>> tfm1.var
        p
        >>> tfm2 = TransferFunctionMatrix(((-S1, -S2), (S1, S2)), (2, 2), p)
        >>> tfm2.var
        p
        >>> tfm3 = TransferFunctionMatrix([[G4]], (1, 1), s)
        >>> tfm3.var
        s

        """
        return self._var

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
        >>> tfm1 = TransferFunctionMatrix((G1, G2), (2, 1), s)
        >>> tfm1.num_inputs
        1
        >>> tfm2 = TransferFunctionMatrix([[G2, -G1, G3], [-G2, -G1, -G3]])
        >>> tfm2.num_inputs
        3

        """
        return self._num_inputs

    @property
    def num_outputs(self):
        """
        Returns the number of outputs of the system.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Series, TransferFunctionMatrix
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> tf3 = TransferFunction(a*s - 4, s**4 + 1, s)
        >>> S1, S2 = Series(tf1, -tf2), Series(tf2, tf3, -tf1)
        >>> tfm1 = TransferFunctionMatrix([tf1, tf2, tf3])
        >>> tfm1.num_outputs
        3
        >>> tfm2 = TransferFunctionMatrix(((S1, S2), (-S2, -S1)))
        >>> tfm2.num_outputs
        2

        """
        return self._num_outputs

    # shape property does not needs to be defined now.

    # @property
    # def shape(self):
    #     """
    #     Returns the shape of the transfer function matrix, that is, ``(# of outputs, # of inputs)``.

    #     If explicitly provided as the second argument of ``TransferFunctionMatrix``, then returns the same,
    #     else, ``TransferFunctionMatrix`` automatically looks up the shape from the list/tuple provided as
    #     the first argument and returns that.

    #     Examples
    #     ========

    #     >>> from sympy.abc import s, p
    #     >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
    #     >>> tf1 = TransferFunction(p**2 - 1, s**4 + s**3 - p, p)
    #     >>> tf2 = TransferFunction(1 - p, p**2 - 3*p + 7, p)
    #     >>> tf3 = TransferFunction(3, 4, p)
    #     >>> tfm1 = TransferFunctionMatrix([[tf1, -tf2]], (1, 2), p)
    #     >>> tfm1.shape
    #     (1, 2)
    #     >>> tfm2 = TransferFunctionMatrix(((-tf2, tf3), (tf1, -tf1)))
    #     >>> tfm2.shape
    #     (2, 2)

    #     """
    #     return self._num_outputs, self._num_inputs

    def __add__(self, other):
        if isinstance(other, (TransferFunctionMatrix, Series)):
            if isinstance(other, Series):
                if other._is_not_matrix:
                    raise ValueError("All the arguments of Series must be either"
                        " TransferFunctionMatrix or Parallel objects.")

            if not self.shape == other.shape:
                raise ShapeError("Shapes of operands are not compatible for addition.")
            if not self.var == other.var:
                raise ValueError("All the TransferFunctionMatrix objects should use the same"
                    " complex variable of the Laplace transform.")
            return Parallel(self, other)

        elif isinstance(other, Parallel):
            if other._is_not_matrix:
                raise ValueError("All the arguments of Parallel must be either TransferFunctionMatrix"
                    " or Series objects.")
            if not self.shape == other.shape:
                raise ShapeError("Shapes of operands are not compatible for addition.")
            if not self.var == other.var:
                raise ValueError("All the TransferFunctionMatrix objects should use the same"
                    " complex variable of the Laplace transform.")
            arg_list = list(other.args)
            return Parallel(self, *arg_list)

        else:
            raise ValueError("TransferFunctionMatrix cannot be added with {}.".
                format(type(other)))

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        # if isinstance(other, (TransferFunctionMatrix, Parallel)):
        #     if isinstance(other, Parallel):
        #         if other._is_not_matrix:
        #             raise ValueError("Only transfer function matrices are allowed in the "
        #                 "parallel configuration.")
        #     if not self.num_inputs == other.num_outputs:
        #         raise ValueError("C = A * B: A has {0} input(s), but B has {1} output(s)."
        #             .format(self.num_inputs, other.num_outputs))
        #     if not self.var == other.var:
        #         raise ValueError("Both TransferFunctionMatrix objects should use the same"
        #             " complex variable of the Laplace transform.")
        #     return Series(self, other)

        # elif isinstance(other, Series):
        #     if other._is_not_matrix:
        #         raise ValueError("Only transfer function matrices are allowed in the "
        #             "series configuration.")
        #     if not self.num_inputs == other.num_outputs:
        #         raise ValueError("C = A * B: A has {0} input(s), but B has {1} output(s)."
        #             .format(self.num_inputs, other.num_outputs))
        #     if not self.var == other.var:
        #         raise ValueError("All the TransferFunctionMatrix objects should use the same"
        #             " complex variable of the Laplace transform.")
        #     arg_list = list(other.args)

        #     return Series(self, *arg_list)
        # else:
        #     raise ValueError("TransferFunctionMatrix cannot be multiplied with {}."
        #         .format(type(other)))

        if isinstance(other, TransferFunctionMatrix):
            if not self.var == other.var:
                raise ValueError("Both TransferFunctionMatrix objects should use the same"
                    " complex variable of the Laplace transform.")
            else:
                # Remember: Series(A, B) = B*A and NOT A*B
                return Series(other, self)
        else:
            raise ValueError("TransferFunctionMatrix cannot be multiplied with {}."
                .format(type(other)))

    def __rmul__(self, other):
        return other*self

    def doit(self, **kwargs):
        """
        Returns the resultant transfer function matrix obtained after evaluating
        the transfer functions in series or parallel configurations (if any present)
        for all entries in a transfer function matrix.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Series, Parallel, TransferFunctionMatrix
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> tf3 = TransferFunction(a*s - 4, s**4 + 1, s)
        >>> S1, S2 = Series(tf1, -tf2), Series(tf2, tf3, -tf1)
        >>> P1, P2 = Parallel(tf2, -tf3), Parallel(-tf1, tf3, -tf2)
        >>> # doit() converts S1, S2, P1, P2 into their resultant transfer functions inside a transfer function matrix.
        >>> tfm1 = TransferFunctionMatrix([S1, S2, -P1])
        >>> tfm1.doit()
        TransferFunctionMatrix([TransferFunction((2 - s**3)*(a*p**2 + b*s), (-p + s)*(s**4 + 5*s + 6), s), TransferFunction((s**3 - 2)*(-a*p**2 - b*s)*(a*s - 4), (-p + s)*(s**4 + 1)*(s**4 + 5*s + 6), s), TransferFunction(-(s**3 - 2)*(s**4 + 1) - (-a*s + 4)*(s**4 + 5*s + 6), (s**4 + 1)*(s**4 + 5*s + 6), s)])
        >>> tfm2 = TransferFunctionMatrix([[tf1, -tf3, -P2]], (1, 3), s)
        >>> tfm2.doit()
        TransferFunctionMatrix([[TransferFunction(a*p**2 + b*s, -p + s, s), TransferFunction(-a*s + 4, s**4 + 1, s), TransferFunction(-(2 - s**3)*(-p + s)*(s**4 + 1) - ((-p + s)*(a*s - 4) + (s**4 + 1)*(-a*p**2 - b*s))*(s**4 + 5*s + 6), (-p + s)*(s**4 + 1)*(s**4 + 5*s + 6), s)]])
        >>> tfm3 = TransferFunctionMatrix(((tf2, P1, P2), (S1, S2, -tf3)), (2, 3), s)
        >>> tfm3.doit()
        TransferFunctionMatrix([[TransferFunction(s**3 - 2, s**4 + 5*s + 6, s), TransferFunction((s**3 - 2)*(s**4 + 1) + (-a*s + 4)*(s**4 + 5*s + 6), (s**4 + 1)*(s**4 + 5*s + 6), s), TransferFunction((2 - s**3)*(-p + s)*(s**4 + 1) + ((-p + s)*(a*s - 4) + (s**4 + 1)*(-a*p**2 - b*s))*(s**4 + 5*s + 6), (-p + s)*(s**4 + 1)*(s**4 + 5*s + 6), s)], [TransferFunction((2 - s**3)*(a*p**2 + b*s), (-p + s)*(s**4 + 5*s + 6), s), TransferFunction((s**3 - 2)*(-a*p**2 - b*s)*(a*s - 4), (-p + s)*(s**4 + 1)*(s**4 + 5*s + 6), s), TransferFunction(-a*s + 4, s**4 + 1, s)]])

        """

        # TODO: Implement doit() for tfm. For now, it returns itself.

        # Function to convert Series(*TransferFunction), Parallel(*TransferFunction) to TransferFunction
        def _to_TransferFunction(arrangement):
            if isinstance(arrangement, (Series, Parallel)):
                return arrangement.doit()
            else:
                return arrangement

        res = self.applyfunc(_to_TransferFunction)
        return res

        # if self.num_inputs == 1:
        #     arg_matrix = list(self.args[0])
        #     for row in range(self.num_outputs):
        #         arg_matrix[row] = arg_matrix[row].doit()
        # else:
        #     arg_matrix = []
        #     for row in range(self.num_outputs):
        #         arg_matrix.append(list(self.args[0][row]))

        #     for row in range(self.num_outputs):
        #         for col in range(self.num_inputs):
        #             arg_matrix[row][col] = arg_matrix[row][col].doit()

        # return TransferFunctionMatrix(arg_matrix)

    # __neg__() is now inherited from ImmutableMatrix class

    # def __neg__(self):
    #     if self.num_inputs == 1:
    #         neg_args = [-col for col in self.args[0]]
    #     else:
    #         neg_args = [[-col for col in row] for row in self.args[0]]
    #     return TransferFunctionMatrix(neg_args)

    # TODO: converts <class 'sympy.physics.control.lti.TransferFunctionMatrix'> to
    # <class 'sympy.matrices.immutable.ImmutableDenseMatrix'> for easing the matrix operations.

    def _to_Immutable_Matrix(self, expr=False):
        self = self.doit()
        _matrix = ImmutableMatrix(self.args[0], self.args[1], list(self.args[2]))
        if not expr:
            return _matrix
        else:
            return _matrix.applyfunc(lambda a: a._to_expr())


    @property
    def is_proper(self):
        """
        Returns True if degree of the numerator polynomial is less than or equal
        to degree of the denominator polynomial for all the SISO transfer functions
        in a transfer function matrix; else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
        >>> tf1 = TransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
        >>> tf2 = TransferFunction(p**2 - 4*p, p**3 + 3*p + 2, p)
        >>> tf3 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf4 = TransferFunction(1, p**2, p)
        >>> tfm1 = TransferFunctionMatrix((-tf1, tf3), (2, 1), s)
        >>> tfm1.is_proper
        False
        >>> tfm2 = TransferFunctionMatrix([[tf2, tf4, -tf2], [-tf4, tf2, tf4]], (2, 3), p)
        >>> tfm2.is_proper
        True

        """
        if self.num_inputs == 1:
            return all(elem.is_proper for elem in self.args[0])
        else:
            return all(self.args[0][row][col].is_proper
                for col in range(self.num_inputs) for row in range(self.num_outputs))

    @property
    def is_strictly_proper(self):
        """
        Returns True if degree of the numerator polynomial is strictly less than
        degree of the denominator polynomial for all the SISO transfer functions
        in a transfer function matrix; else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> tf3 = TransferFunction(a, s**2 + b, s)
        >>> tfm1 = TransferFunctionMatrix([-tf1, tf2, -tf3])
        >>> tfm1.is_strictly_proper
        False
        >>> tfm2 = TransferFunctionMatrix([[tf2, tf3], [-tf3, -tf2]])
        >>> tfm2.is_strictly_proper
        True

        """
        if self.num_inputs == 1:
            return all(elem.is_strictly_proper for elem in self.args[0])
        else:
            return all(self.args[0][row][col].is_strictly_proper
                for col in range(self.num_inputs) for row in range(self.num_outputs))

    @property
    def is_biproper(self):
        """
        Returns True if degree of the numerator polynomial is equal to
        degree of the denominator polynomial for all the SISO transfer
        functions in a transfer function matrix; else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(-b*p**4 - a*s**2, a - s**2, s)
        >>> tf3 = TransferFunction(s**2, s + a, s)
        >>> tfm1 = TransferFunctionMatrix([tf1, -tf2])
        >>> tfm1.is_biproper
        True
        >>> tfm2 = TransferFunctionMatrix(((tf1, tf2), (tf3, -tf2)))
        >>> tfm2.is_biproper
        False

        """
        if self.num_inputs == 1:
            return all(elem.is_biproper for elem in self.args[0])
        else:
            return all(self.args[0][row][col].is_biproper
                for col in range(self.num_inputs) for row in range(self.num_outputs))
