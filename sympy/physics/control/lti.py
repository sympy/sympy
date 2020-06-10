from sympy import Basic, Mul, degree, Symbol, expand, cancel, Expr
from sympy.core.numbers import Integer
from sympy.core.sympify import sympify, _sympify

__all__ = ['TransferFunction',]


class TransferFunction(Basic):
    """TransferFunction(num, den, var)

    A class for representing Transfer Functions. This class is used to represent
    LTI (Linear, time-invariant) systems in transfer function form. The arguments
    are ``num``, ``den``, and ``var``, where ``num`` and ``den`` are numerator and
    denominator of the ``TransferFunction`` respectively, and the third argument is
    a variable used to anchor the symbols. ``num`` and ``den`` can be either
    sympy expressions or numbers, whereas ``var`` has to be a Symbol. This class is a
    candidate for a full symbolic model rather than just being a simple solver.

    Parameters
    ==========

    num : Numerator
        Numerator of the Transfer Function.
    den : Denominator
        Denominator of the Transfer Function.
    var : Variable
        Variable used to anchor the symbols.

    Raises
    ======

    TypeError
        When `var` is not a Symbol or when incorrect type for `num` or `den`
        is supplied.
    ValueError
        When `den` is zero.

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

    Using different polynomial variables.

    >>> tf2 = TransferFunction(a*p**3 - a*p**2 + s*p, p + a**2, p)
    >>> tf2
    TransferFunction(a*p**3 - a*p**2 + p*s, a**2 + p, p)
    >>> tf3 = TransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p)
    >>> tf3
    TransferFunction((p - 1)*(p + 3), (p - 1)*(p + 5), p)

    Using Float or Integer (or other constants) as numerator and denominator.

    >>> tf4 = TransferFunction(1/2, 4, s)
    >>> tf4.num
    0.500000000000000
    >>> tf4.den
    4
    >>> tf4.var
    s

    References
    ==========

    Joao P. Hespanha, Linear Systems Theory. 2009.

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
        """Returns the Numerator of the Transfer Function.

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
        """Returns the Denominator of the Transfer Function.

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
        """Returns the variable that is used to anchor the symbols.

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

    @property
    def bound_symbols(self):
        """Returns only variables that are dummy variables.

        Examples
        ========

        >>> from sympy.abc import s, p, l
        >>> from sympy.physics.control.lti import TransferFunction
        >>> G1 = TransferFunction(s + l, s - 2*p, s)
        >>> G1.bound_symbols
        [s]
        >>> G2 = TransferFunction(p**2 + 2*p + 4, p - 6, p)
        >>> G2.bound_symbols
        [p]

        """
        return [self.var]

    def _eval_subs(self, old, new):
        arg_num = self.args[0].subs(old, new)
        arg_den = self.args[1].subs(old, new)
        argnew = TransferFunction(arg_num, arg_den, self.var)
        return self if old in self.bound_symbols else argnew

    def _eval_simplify(self, **kwargs):
        tf = cancel(Mul(self.num, 1/self.den)).as_numer_denom()
        num_, den_ = tf[0], tf[1]
        return TransferFunction(num_, den_, self.var)

    def expand(self):
        """
        Returns the Transfer Function with Numerator and Denominator
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
        """Adds Transfer Functions.

        Parameters
        ==========

        other : TransferFunction
            The Transfer Function to add to current (self) Transfer Function.

        Returns
        =======

        TransferFunction
            The resultant Transfer Function after adding (with '+' operator) self to other.

        Examples
        ========

        >>> from sympy.abc import s, a
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction(s + 6, s - 5, s)
        >>> tf2 = TransferFunction(s + 3, s + 1, s)
        >>> tf1 + tf2
        TransferFunction((s - 5)*(s + 3) + (s + 1)*(s + 6), (s - 5)*(s + 1), s)
        >>> tf1 + (a - 1)
        TransferFunction(s + (a - 1)*(s - 5) + 6, s - 5, s)
        >>> tf1 + 8
        TransferFunction(9*s - 34, s - 5, s)

        """
        other = _sympify(other)
        if other.is_number:
            return TransferFunction(self.num + self.den*other, self.den, self.var)

        elif isinstance(other, TransferFunction):
            if not self.var == other.var:
                raise ValueError("Both the Transfer Functions should be anchored "
                    "with the same variable.")
            p = self.num * other.den + other.num * self.den
            q = self.den * other.den
            return TransferFunction(p, q, self.var)

        elif isinstance(other, Expr) and other.has(Symbol):
            # other input is a polynomial.
            if not other.is_commutative:
                raise ValueError("Only commutative expressions can be added "
                    "with a TransferFunction.")
            return TransferFunction(self.num + self.den*other, self.den, self.var)
        else:
            raise ValueError("TransferFunction cannot be added with {}.".
                format(type(other)))

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        """Subtracts Transfer Functions.

        Parameters
        ==========

        other : TransferFunction
            The Transfer Function to subtract from current (self) Transfer Function.

        Returns
        =======

        TransferFunction
            The resultant Transfer Function after subtracting (with '-' operator)
            other from self.

        Examples
        ========

        >>> from sympy.abc import p, b
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction(p**2, p - 4, p)
        >>> tf2 = TransferFunction(p + 1, 7 - p, p)
        >>> tf1 - tf2
        TransferFunction(p**2*(7 - p) - (p - 4)*(p + 1), (7 - p)*(p - 4), p)
        >>> tf1 - 4
        TransferFunction(p**2 - 4*p + 16, p - 4, p)
        >>> tf2 - (b + 7)
        TransferFunction(p - (7 - p)*(b + 7) + 1, 7 - p, p)

        """
        other = _sympify(other)
        if other.is_number:
            return TransferFunction(self.num - self.den*other, self.den, self.var)

        elif isinstance(other, TransferFunction):
            if not self.var == other.var:
                raise ValueError("Both the Transfer Functions should be anchored "
                    "with the same variable.")
            p = self.num * other.den - other.num * self.den
            q = self.den * other.den
            return TransferFunction(p, q, self.var)

        elif isinstance(other, Expr) and other.has(Symbol):
            return TransferFunction(self.num - self.den*other, self.den, self.var)
        else:
            raise ValueError("{} cannot be subtracted from TransferFunction."
                .format(type(other)))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        """Multiplies Transfer Functions.

        Parameters
        ==========

        other : TransferFunction
            The Transfer Function to multiply to current (self) Transfer Function.

        Returns
        =======

        TransferFunction
            The resultant Transfer Function after multiplying (with '*' operator)
            self to other.

        Examples
        ========

        >>> from sympy.abc import s, p, a
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction(s + 1, s - 8, s)
        >>> tf2 = TransferFunction(p + s**2, p - s, s)
        >>> tf1*tf2
        TransferFunction((p + s**2)*(s + 1), (p - s)*(s - 8), s)
        >>> 9*tf2
        TransferFunction(9*p + 9*s**2, p - s, s)
        >>> tf1*(a + 6)
        TransferFunction((a + 6)*(s + 1), s - 8, s)

        """
        other = _sympify(other)
        if other.is_number:
            return TransferFunction(self.num*other, self.den, self.var)

        elif isinstance(other, TransferFunction):
            if not self.var == other.var:
                raise ValueError("Both the Transfer Functions should be anchored "
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

    def __div__(self, other):
        """Divides Transfer Functions.

        Parameters
        ==========

        other : TransferFunction
            The Transfer Function to be divided from current (self) Transfer Function.

        Returns
        =======

        TransferFunction
            Resultant Transfer Function after dividing (with '/' operator) self by other.

        Examples
        ========

        >>> from sympy.abc import s, p, a
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction(s + 1, s - 8, s)
        >>> tf2 = TransferFunction(p + s**2, p - s, s)
        >>> tf1/tf2
        TransferFunction((p - s)*(s + 1), (p + s**2)*(s - 8), s)
        >>> tf2/(-3)
        TransferFunction(p + s**2, -3*p + 3*s, s)
        >>> tf1/(a**2)
        TransferFunction(s + 1, a**2*(s - 8), s)

        """
        other = _sympify(other)
        if other.is_number:
            if other == 0:
                raise ValueError("TransferFunction cannot be divided by zero.")
            return TransferFunction(self.num, self.den*other, self.var)

        elif isinstance(other, TransferFunction):
            if not self.var == other.var:
                raise ValueError("Both the Transfer Functions should be anchored "
                    "with the same variable.")
            p = self.num * other.den
            q = self.den * other.num
            return TransferFunction(p, q, self.var)

        elif isinstance(other, Expr) and other.has(Symbol):
            return TransferFunction(self.num, self.den*other, self.var)
        else:
            raise ValueError("TransferFunction cannot be divided by {}.".
                format(type(other)))

    __truediv__ = __div__

    def __rtruediv__(self, other):
        return _sympify(other) * self**-1

    __rdiv__ = __rtruediv__

    def __pow__(self, p):
        """Finds the p-th power of the Transfer Function.

        Parameters
        ==========

        p : int
            Power to be applied on the Transfer Function.

        Returns
        =======

        TransferFunction
            Returns the p-th power of the current Transfer Function.

        Examples
        ========

        >>> from sympy.abc import s, p, a
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction(s + a, s - a, s)
        >>> tf1**3
        TransferFunction((a + s)**3, (-a + s)**3, s)
        >>> tf2 = TransferFunction(s, s**2 + p, s)
        >>> tf2**0
        TransferFunction(1, 1, s)
        >>> tf3 = TransferFunction(p + 4, p - 3, p)
        >>> tf3**-1
        TransferFunction(p - 3, p + 4, p)

        """
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
        """
        Negates a Transfer Function.

        Examples
        ========

        >>> from sympy.abc import s, p, a
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction(-a + s, p**2 + s, p)
        >>> -tf1
        TransferFunction(a - s, p**2 + s, p)
        >>> tf2 = TransferFunction(s**4 - 2*s**3 + 5*s + 4, s + 4, s)
        >>> -tf2
        TransferFunction(-s**4 + 2*s**3 - 5*s - 4, s + 4, s)

        """
        return TransferFunction(-self.num, self.den, self.var)

    @property
    def is_proper(self):
        """
        Returns True if degree of the Numerator is less than or equal to
        degree of the Denominator, else False.

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
        Returns True if degree of the Numerator is strictly less than
        degree of the Denominator, else False.

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
        Returns True if degree of the Numerator is equal to degree
        of the Denominator, else False.

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
