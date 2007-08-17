"""Class for modular integer arithmetic."""

class ModularInteger(object):
    """Modular integer arithmetic, based on Python integers."""

    modulus = 0

    def __init__(self, value):
        self.value = value % self.modulus

    def __repr__(self):
        return "%s(%s)" % (self.__class__, self.value)

    def __str__(self):
        return "%s mod %s" % (int(self), self.modulus)

    def __int__(self):
        """Return the unique integer -m/2 < i <= m/2."""
        if self.value <= self.modulus/2:
            return self.value
        else:
            return self.value - self.modulus

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.value)

    def __add__(self, other):
        return self.__class__(self.value + other.value)

    def __sub__(self, other):
        return self.__class__(self.value - other.value)

    def __mul__(self, other):
        return self.__class__(self.value * other.value)

    def __div__(self, other):
        g, x, y = xgcd(other.value, self.modulus)
        assert g == 1, "Zero division!"
        return self.__class__(self.value * x)

    def __pow__(self, exponent):
        """Repeated squaring."""
        exponent = int(exponent)
        if exponent < 0:
            g, x, y = xgcd(self.value, self.exponent)
            assert g == 1, "Zero division!"
            value = x
            exponent *= -1
        if exponent == 0:
            return self.__class__(1)
        else:
            value = self.value
        value = pow(value, exponent, self.modulus)
        return self.__class__(value)

    def __eq__(self, other):
        return self.value == other.value

    def __ne__(self, other):
        return self.value != other.value

    def __nonzero__(self):
        return bool(self.value)

def ModularIntegerFactory(m):
    """Create custom class for specific integer modulus."""

    class newClass(ModularInteger):
        modulus = m
        
    newClass.__name__ = "IntMod%s" % m
    return newClass

def gcd(a, b):
    if a == 0: return abs(b)
    if b == 0: return abs(a)
    if b < 0: b *= -1
    while b:
        a, b = b, a % b
    return a

def xgcd(a, b):
    """
    Returns g, x, y such that g = x*a + y*b = gcd(a,b).
    Input:
        a -- an integer
        b -- an integer
    Output:
        g -- an integer, the gcd of a and b
        x -- an integer
        y -- an integer
    Examples:
    >>> xgcd(2,3)
    (1, -1, 1)
    >>> xgcd(10, 12)
    (2, -1, 1)
    >>> g, x, y = xgcd(100, 2004)
    >>> print g, x, y
    4 -20 1
    >>> print x*100 + y*2004
    4
    """
    if a == 0 and b == 0: return (0, 0, 1)
    if a == 0: return (abs(b), 0, b/abs(b))
    if b == 0: return (abs(a), a/abs(a), 0)
    x_sign = 1; y_sign = 1
    if a < 0: a = -a; x_sign = -1
    if b < 0: b = -b; y_sign = -1
    x = 1; y = 0; r = 0; s = 1
    while b != 0:
        (c, q) = (a%b, a/b)
        (a, b, r, s, x, y) = (b, c, x-q*r, y-q*s, r, s)
    return (a, x*x_sign, y*y_sign)

def crt(m, v, symmetric=False):
    """Chinese remainder theorem.

    The integers in m are assumed to be pairwise coprime. The output
    is then an integer f, such that that f = v_i mod m_i for each pair
    out of v and m.
    """
    mm = 1
    for m_i in m:
        mm *= m_i
    result = 0
    for m_i, v_i in zip(m, v):
        e = mm/m_i
        g, s, t = xgcd(e, m_i)
        c = (v_i*s) % m_i
        result += c*e
    result %= mm
    if symmetric:
        if result <= mm/2:
            return result
        else:
            return result - mm
    else:
        return result

def crt1(m):
    """First part of chines remainder theorem, for multiple application."""
    mm = 1
    e = []
    s = []
    for m_i in m:
        mm *= m_i
    for m_i in m:
        e.append(mm/m_i)
        s.append(xgcd(e[-1], m_i)[1])
    return mm, e, s

def crt2(m, v, mm, e, s, symmetric=False):
    """Second part of chines remainder theorem, for multiple application."""
    result = 0
    for m_i, v_i, e_i, s_i in zip(m, v, e, s):
        c = v_i*s_i % m_i
        result += c*e_i
    result %= mm
    if symmetric:
        if result <= mm/2:
            return result
        else:
            return result - mm
    else:
        return result
    
