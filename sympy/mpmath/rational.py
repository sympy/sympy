# TODO: use gmpy.mpq when available?

class mpq(tuple):
    """
    Rational number type, only intended for internal use.
    """

    """
    def _mpmath_(self, prec, rounding):
        # XXX
        return mp.make_mpf(from_rational(self[0], self[1], prec, rounding))
        #(mpf(self[0])/self[1])._mpf_

    """

    def __int__(self):
        a, b = self
        return a // b

    def __abs__(self):
        a, b = self
        return mpq((abs(a), b))

    def __neg__(self):
        a, b = self
        return mpq((-a, b))

    def __nonzero__(self):
        return bool(self[0])

    def __cmp__(self, other):
        if type(other) is int and self[1] == 1:
            return cmp(self[0], other)
        return NotImplemented

    def __add__(self, other):
        if isinstance(other, mpq):
            a, b = self
            c, d = other
            return mpq((a*d+b*c, b*d))
        if isinstance(other, (int, long)):
            a, b = self
            return mpq((a+b*other, b))
        return NotImplemented

    __radd__ = __add__

    def __sub__(self, other):
        if isinstance(other, mpq):
            a, b = self
            c, d = other
            return mpq((a*d-b*c, b*d))
        if isinstance(other, (int, long)):
            a, b = self
            return mpq((a-b*other, b))
        return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, mpq):
            a, b = self
            c, d = other
            return mpq((b*c-a*d, b*d))
        if isinstance(other, (int, long)):
            a, b = self
            return mpq((b*other-a, b))
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, mpq):
            a, b = self
            c, d = other
            return mpq((a*c, b*d))
        if isinstance(other, (int, long)):
            a, b = self
            return mpq((a*other, b))
        return NotImplemented

    def __div__(self, other):
        if isinstance(other, (int, long)):
            if other:
                a, b = self
                return mpq((a, b*other))
            raise ZeroDivisionError
        return NotImplemented

    def __pow__(self, other):
        if type(other) is int:
            a, b = self
            return mpq((a**other, b**other))
        return NotImplemented

    __rmul__ = __mul__


mpq_1 = mpq((1,1))
mpq_0 = mpq((0,1))
mpq_1_2 = mpq((1,2))
mpq_3_2 = mpq((3,2))
mpq_1_4 = mpq((1,4))
mpq_1_16 = mpq((1,16))
mpq_3_16 = mpq((3,16))
mpq_5_2 = mpq((5,2))
mpq_3_4 = mpq((3,4))
mpq_7_4 = mpq((7,4))
mpq_5_4 = mpq((5,4))

