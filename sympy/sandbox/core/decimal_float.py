
import math
import decimal_math
from utils import memoizer_immutable_args
from basic import Basic, sympify
from number import Real

import decimal
from decimal import Decimal

class Float(Real, Decimal):

    @memoizer_immutable_args('Float.__new__')
    def __new__(cls, f):
        if isinstance(f, Basic):
            return f.evalf()
        if isinstance(f, float):
            #f = str(f)
            f = cls.float_to_decimal(f)
        return Decimal.__new__(cls, f)

    def __float__(self):
        return float(self.as_native())

    def as_native(self):
        return Decimal(self)

    @staticmethod
    def float_to_decimal(f):
        "Convert a floating point number to a Decimal with no loss of information"
        # Transform (exactly) a float to a mantissa (0.5 <= abs(m) < 1.0) and an
        # exponent.  Double the mantissa until it is an integer.  Use the integer
        # mantissa and exponent to compute an equivalent Decimal.  If this cannot
        # be done exactly, then retry with more precision.

        try:
            mantissa, exponent = math.frexp(f)
        except OverflowError:
            return decimal.Inf

        while mantissa != int(mantissa):
            mantissa *= 2.0
            exponent -= 1
        mantissa = int(mantissa)

        oldcontext = decimal.getcontext()
        decimal.setcontext(decimal.Context(traps=[decimal.Inexact]))
        try:
            while True:
                try:
                    return mantissa * Decimal(2) ** exponent
                except decimal.Inexact:
                    decimal.getcontext().prec += 1
        finally:
            decimal.setcontext(oldcontext)

    def __pow__(self, other):
        other = sympify(other).evalf()
        if other.is_Float:
            r = self._eval_power(other)
            if r is not None:
                return r
        return Basic.Pow(self, other)

    def _eval_power(b, e):
        """
        b is Real but not equal to rationals, integers, 0.5, oo, -oo, nan
        e is symbolic object but not equal to 0, 1

        (-p) ** r -> exp(r * log(-p)) -> exp(r * (log(p) + I*Pi)) ->
                  -> p ** r * (sin(Pi*r) + cos(Pi*r) * I)
        """
        if e.is_Number:
            b2 = b.as_native()
            if e.is_Integer:
                e = int(e.p)
                return Real(decimal_math.pow(b2, e))
            e = e.evalf()
            e2 = e.as_native()
            if b.is_negative and not e.is_integer:
                m = decimal_math.pow(-b2, e2)
                a = decimal_math.pi() * e2
                s = m * decimal_math.sin(a)
                c = m * decimal_math.cos(a)
                return Float(s) + Float(c) * Basic.ImaginaryUnit()
            return Float(decimal_math.pow(b2, e2))
        return
