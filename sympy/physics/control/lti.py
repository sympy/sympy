from sympy import Basic, Mul, degree, Symbol, expand, cancel, Expr, exp
from sympy.core.evalf import EvalfMixin
from sympy.core.numbers import Integer
from sympy.core.sympify import sympify, _sympify

__all__ = ['TransferFunction', 'Series', 'Parallel', 'Feedback']


class TransferFunction(Basic, EvalfMixin):
    """TransferFunction(num, den, var)

    A class for representing continuous transfer functions. This class is used to represent
    LTI (Linear, time-invariant) systems in transfer function form. The arguments
    are ``num``, ``den``, and ``var``, where ``num`` and ``den`` are numerator and
    denominator polynomials of the ``TransferFunction`` respectively, and the third argument is
    a variable used to anchor these polynomials of the transfer function. ``num`` and ``den`` can
    be either sympy expressions (with no time delay terms) or numbers, whereas ``var``
    has to be a Symbol.

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
        a number or expression/polynomial. Also, when ``num`` or ``den`` has
        a time delay term.
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
    >>> tf7**0
    TransferFunction(1, 1, s)
    >>> tf8 = TransferFunction(p + 4, p - 3, p)
    >>> tf8**-1
    TransferFunction(p - 3, p + 4, p)

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
        if isinstance(other, (TransferFunction, Series)):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
            return Parallel(self, other)
        elif isinstance(other, Parallel):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
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
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
            return Parallel(self, -other)
        elif isinstance(other, Parallel):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
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
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
            return Series(self, other)
        elif isinstance(other, Series):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
            arg_list = list(other.args)
            return Series(self, *arg_list)
        else:
            raise ValueError("TransferFunction cannot be multiplied with {}."
                .format(type(other)))

    __rmul__ = __mul__

    def __truediv__(self, other):
        if (isinstance(other, Parallel) and isinstance(other.args[0], TransferFunction)
            and isinstance(other.args[1], (Series, TransferFunction))):

            if not self.var == other.var:
                raise ValueError("Both TransferFunction and Parallel should be anchored "
                    "with the same variable.")
            if other.args[1] == self:
                return Feedback(self, other.args[0])
            other_arg_list = list(other.args[1].args)
            if self in other_arg_list:
                other_arg_list.remove(self)
            else:
                return Feedback(self, Series(*other_arg_list))

            if len(other_arg_list) == 0:
                return Feedback(self, other.args[0])
            elif len(other_arg_list) == 1:
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


class Series(Basic):

    def __new__(cls, *args, evaluate=False):
        if not all(isinstance(arg, (TransferFunction, Parallel, Series)) for arg in args):
            raise TypeError("Unsupported type of argument(s) for Series.")

        obj = super(Series, cls).__new__(cls, *args)
        obj._var = None
        for arg in args:
            if obj._var is None:
                obj._var = arg.var
            elif obj._var != arg.var:
                raise ValueError("All transfer functions should be anchored with the same variable.")
        if evaluate:
            return obj.doit()
        return obj

    @property
    def var(self):
        return self._var

    def doit(self, **kwargs):
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
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")

            return Parallel(self, other)
        elif isinstance(other, Parallel):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
            arg_list = list(other.args)

            return Parallel(self, *arg_list)
        else:
            raise ValueError("This transfer function expression is invalid.")

    def __sub__(self, other):
        if isinstance(other, (TransferFunction, Series)):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")

            return Parallel(self, -other)
        elif isinstance(other, Parallel):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
            arg_list = [-i for i in list(other.args)]

            return Parallel(self, *arg_list)
        else:
            raise ValueError("This transfer function expression is invalid.")

    def __mul__(self, other):
        if isinstance(other, (TransferFunction, Parallel)):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
            arg_list = list(self.args)

            return Series(*arg_list, other)
        elif isinstance(other, Series):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
            self_arg_list = list(self.args)
            other_arg_list = list(other.args)

            return Series(*self_arg_list, *other_arg_list)
        else:
            raise ValueError("This transfer function expression is invalid.")

    def __truediv__(self, other):
        if (isinstance(other, Parallel) and isinstance(other.args[0], TransferFunction)
            and isinstance(other.args[1], Series)):

            if not self.var == other.var:
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
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
        return self.doit().is_proper

    @property
    def is_strictly_proper(self):
        return self.doit().is_strictly_proper

    @property
    def is_biproper(self):
        return self.doit().is_biproper


class Parallel(Basic):

    def __new__(cls, *args, evaluate=False):
        if not all(isinstance(arg, (TransferFunction, Series, Parallel)) for arg in args):
            raise TypeError("Unsupported type of argument(s) for Parallel.")

        obj = super(Parallel, cls).__new__(cls, *args)
        obj._var = None
        for arg in args:
            if obj._var is None:
                obj._var = arg.var
            elif obj._var != arg.var:
                raise ValueError("All transfer functions should be anchored with the same variable.")
        if evaluate:
            return obj.doit()
        return obj

    @property
    def var(self):
        return self._var

    def doit(self, **kwargs):
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
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
            arg_list = list(self.args)
            arg_list.append(other)

            return Parallel(*arg_list)
        elif isinstance(other, Parallel):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
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
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
            arg_list = list(self.args)
            arg_list.append(-other)

            return Parallel(*arg_list)
        elif isinstance(other, Parallel):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
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
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
            return Series(self, other)
        elif isinstance(other, Series):
            if not self.var == other.var:
                raise ValueError("All the transfer functions should be anchored "
                    "with the same variable.")
            arg_list = list(other.args)
            return Series(self, *arg_list)
        else:
            raise ValueError("This transfer function expression is invalid.")

    def __neg__(self):
        return Series(TransferFunction(-1, 1, self.var), self)

    @property
    def is_proper(self):
        return self.doit().is_proper

    @property
    def is_strictly_proper(self):
        return self.doit().is_strictly_proper

    @property
    def is_biproper(self):
        return self.doit().is_biproper


class Feedback(Basic):

    def __new__(cls, num, den):
        if not (isinstance(num, (TransferFunction, Series)) and num != den
            and isinstance(den, (TransferFunction, Series))):
            raise TypeError("Unsupported type for numerator or denominator of Feedback.")

        if not num.var == den.var:
            raise ValueError("Both numerator and denominator should be anchored with the same variable.")
        obj = super(Feedback, cls).__new__(cls, num, den)
        obj._num = num
        obj._den = den
        obj._var = num.var

        return obj

    @property
    def num(self):
        return self._num

    @property
    def den(self):
        return self._den

    @property
    def var(self):
        return self._var

    def doit(self, **kwargs):
        arg_list = list(self.num.args) if isinstance(self.num, Series) else [self.num]
        # F_n and F_d are resultant TFs of num and den of Feedback.
        F_n, tf = self.num.doit(), TransferFunction(1, 1, self.num.var)
        F_d = Parallel(tf, Series(self.den, *arg_list)).doit()

        return TransferFunction(F_n.num*F_d.den, F_n.den*F_d.num, F_n.var)

    def _eval_rewrite_as_TransferFunction(self, num, den, **kwargs):
        return self.doit()

    def __neg__(self):
        return Feedback(-self.num, self.den)
