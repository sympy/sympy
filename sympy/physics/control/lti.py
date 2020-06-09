from sympy import Basic, Mul, degree, Symbol, expand, cancel, Expr
from sympy.core.numbers import Integer, Float
from sympy.core.sympify import sympify, _sympify

__all__ = ['TransferFunction',]


class TransferFunction(Basic):

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
            raise ValueError("Unsupported type for numerator or denominator of TransferFunction.")

    @property
    def num(self):
        return self._num

    @property
    def den(self):
        return self._den

    @property
    def var(self):
        return self._var

    @property
    def bound_symbols(self):
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
        return TransferFunction(expand(self.num), expand(self.den), self.var)

    def __add__(self, other):
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
        return degree(self.num) <= degree(self.den)

    @property
    def is_strictly_proper(self):
        return degree(self.num) < degree(self.den)

    @property
    def is_biproper(self):
        return degree(self.num) == degree(self.den)
