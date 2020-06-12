from sympy import Basic, Mul, degree, Symbol, expand, cancel, Expr
from sympy.core.evalf import EvalfMixin
from sympy.core.numbers import Integer
from sympy.core.sympify import sympify, _sympify

__all__ = ['TransferFunction',]


class TransferFunction(Basic, EvalfMixin):
    """TransferFunction(num, den, var)

    A class for representing continuous transfer functions. This class is used to represent
    LTI (Linear, time-invariant) systems in transfer function form. The arguments
    are ``num``, ``den``, and ``var``, where ``num`` and ``den`` are numerator and
    denominator polynomials of the ``TransferFunction`` respectively, and the third argument is
    a variable used to anchor these polynomials of the transfer function. ``num`` and ``den`` can
    be either sympy expressions or numbers, whereas ``var`` has to be a Symbol.

    Parameters
    ==========

    num : Expr, Number
        The numerator polynomial of the transfer function.
    den : Expr, Number
        The denominator polynomial of the transfer function.
    var : Symbol
        Variable used to anchor the polynomials of the transfer function.

    Raises
    ======

    TypeError
        When ``var`` is not a Symbol or when ``num`` or ``den`` is not
        a number or expression/polynomial.
    ValueError
        When ``den`` is zero.

    Examples
    ========

    >>> from sympy.abc import s, p, a, b
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

    Using different polynomial variables.

    >>> tf2 = TransferFunction(a*p**3 - a*p**2 + s*p, p + a**2, p)
    >>> tf2
    TransferFunction(a*p**3 - a*p**2 + p*s, a**2 + p, p)
    >>> tf3 = TransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p)
    >>> tf3
    TransferFunction((p - 1)*(p + 3), (p - 1)*(p + 5), p)

    Using ``-`` operator in front of a transfer function to negate it.

    >>> tf4 = TransferFunction(-a + s, p**2 + s, p)
    >>> -tf4
    TransferFunction(a - s, p**2 + s, p)
    >>> tf5 = TransferFunction(s**4 - 2*s**3 + 5*s + 4, s + 4, s)
    >>> -tf5
    TransferFunction(-s**4 + 2*s**3 - 5*s - 4, s + 4, s)

    Using Float or Integer (or other constants) as numerator and denominator.

    >>> tf6 = TransferFunction(1/2, 4, s)
    >>> tf6.num
    0.500000000000000
    >>> tf6.den
    4
    >>> tf6.var
    s

    Using ``**`` operator to take the power of a transfer function.

    >>> tf7 = TransferFunction(s + a, s - a, s)
    >>> tf7**3
    TransferFunction((a + s)**3, (-a + s)**3, s)
    >>> tf7 = TransferFunction(s, s**2 + p, s)
    >>> tf8**0
    TransferFunction(1, 1, s)
    >>> tf8 = TransferFunction(p + 4, p - 3, p)
    >>> tf8**-1
    TransferFunction(p - 3, p + 4, p)


    We can add the current (self) transfer function to a polynomial, a number or another
    transfer function using ``+`` operator. This can be called as the parallel inter-connection of
    two transfer functions.

    >>> G1 = TransferFunction(s + 6, s - 5, s)
    >>> G2 = TransferFunction(s + 3, s + 1, s)
    >>> G1 + G2
    TransferFunction((s - 5)*(s + 3) + (s + 1)*(s + 6), (s - 5)*(s + 1), s)
    >>> G1 + (a - 1)
    TransferFunction(s + (a - 1)*(s - 5) + 6, s - 5, s)
    >>> G1 + 8
    TransferFunction(9*s - 34, s - 5, s)

    We can subtract a polynomial, a number or a transfer function from the current (self)
    transfer function using ``-`` operator.

    >>> G1 - G2
    TransferFunction((s + 1)*(s + 6) - (s - 5)*(s + 3), (s - 5)*(s + 1), s)
    >>> G1 - 4
    TransferFunction(26 - 3*s, s - 5, s)
    >>> G2 - (b + 7)
    TransferFunction(s - (b + 7)*(s + 1) + 3, s + 1, s)

    We can also multiply the current (self) transfer function with a polynomial, a number or
    another transfer function using ``*`` operator. This can be called as the series inter-connection
    of two transfer functions.

    >>> G3 = TransferFunction(s + 1, s - 8, s)
    >>> G4 = TransferFunction(p + s**2, p - s, s)
    >>> G3*G4
    TransferFunction((p + s**2)*(s + 1), (p - s)*(s - 8), s)
    >>> 9*G4
    TransferFunction(9*p + 9*s**2, p - s, s)
    >>> G3*(a + 6)
    TransferFunction((a + 6)*(s + 1), s - 8, s)

    Similarly, we can divide the current (self) transfer function by a polynomial, a number or
    another transfer function using ``/`` operator.

    >>> G3/G4
    TransferFunction((p - s)*(s + 1), (p + s**2)*(s - 8), s)
    >>> G4/(-3)
    TransferFunction(p + s**2, -3*p + 3*s, s)
    >>> G3/(a**2)
    TransferFunction(s + 1, a**2*(s - 8), s)

    References
    ==========

    .. [1] http://www.cds.caltech.edu/~murray/amwiki/index.php?title=First_Edition

    """
    def __new__(cls, num, den, var):
        num, den = _sympify(num), _sympify(den)

        if not isinstance(var, Symbol):
            raise TypeError("Variable input must be a Symbol.")
        if den == 0:
            raise ValueError("TransferFunction can't have a zero denominator.")

        if (((isinstance(num, Expr) and num.has(Symbol)) or num.is_number) and
            ((isinstance(den, Expr) and den.has(Symbol)) or den.is_number)):
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
        Returns the numerator polynomial/expression of the transfer function.

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
        Returns the denominator polynomial/expression of the transfer function.

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
        Returns the variable used to anchor the expressions/polynomials of
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

    def __add__(self, other):
        other = _sympify(other)
        if other.is_number:
            return TransferFunction(self.num + self.den*other, self.den, self.var)

        elif isinstance(other, TransferFunction):
            if not self.var == other.var:
                raise ValueError("Both the transfer functions should be anchored "
                    "with the same variable.")
            p = self.num * other.den + other.num * self.den
            q = self.den * other.den
            return TransferFunction(p, q, self.var)

        elif isinstance(other, Expr) and other.has(Symbol):
            # other input is a polynomial.
            if not other.is_commutative:
                raise ValueError("Only commutative expressions can be added "
                    "with a transfer function.")
            return TransferFunction(self.num + self.den*other, self.den, self.var)
        else:
            raise ValueError("TransferFunction cannot be added with {}.".
                format(type(other)))

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        other = _sympify(other)
        if other.is_number:
            return TransferFunction(self.num - self.den*other, self.den, self.var)

        elif isinstance(other, TransferFunction):
            if not self.var == other.var:
                raise ValueError("Both the transfer functions should be anchored "
                    "with the same variable.")
            p = self.num * other.den - other.num * self.den
            q = self.den * other.den
            return TransferFunction(p, q, self.var)

        elif isinstance(other, Expr) and other.has(Symbol):
            return TransferFunction(self.num - self.den*other, self.den, self.var)
        else:
            raise ValueError("{} cannot be subtracted from a TransferFunction."
                .format(type(other)))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        other = _sympify(other)
        if other.is_number:
            return TransferFunction(self.num*other, self.den, self.var)

        elif isinstance(other, TransferFunction):
            if not self.var == other.var:
                raise ValueError("Both the transfer functions should be anchored "
                    "with the same variable.")
            p = self.num * other.num
            q = self.den * other.den
            return TransferFunction(p, q, self.var)

        elif isinstance(other, Expr) and other.has(Symbol):
            # other input is a polynomial.
            if not other.is_commutative:
                raise ValueError("Only commutative expressions can be multiplied "
                    "with a TransferFunction.")
            return TransferFunction(self.num*other, self.den, self.var)
        else:
            raise ValueError("TransferFunction cannot be multiplied with {}."
                .format(type(other)))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = _sympify(other)
        if other.is_number:
            if other == 0:
                raise ValueError("TransferFunction cannot be divided by zero.")
            return TransferFunction(self.num, self.den*other, self.var)

        elif isinstance(other, TransferFunction):
            if not self.var == other.var:
                raise ValueError("Both the transfer functions should be anchored "
                    "with the same variable.")
            p = self.num * other.den
            q = self.den * other.num
            return TransferFunction(p, q, self.var)

        elif isinstance(other, Expr) and other.has(Symbol):
            return TransferFunction(self.num, self.den*other, self.var)
        else:
            raise ValueError("TransferFunction cannot be divided by {}.".
                format(type(other)))

    def __rtruediv__(self, other):
        return _sympify(other) * self**-1

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
