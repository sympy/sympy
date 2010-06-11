"""Classes defining properties of ground domains, e.g. ZZ, QQ, ZZ[x] ... """

from sympy.core import S, Basic, sympify

from sympy.polys.polyerrors import (
    ExactQuotientFailed,
    IsomorphismFailed,
    UnificationFailed,
    GeneratorsNeeded,
    CoercionFailed,
    NotInvertible,
    NotAlgebraic,
    DomainError,
)

import math

def factorial(m):
    """
    Returns `m!`.

    Example
    -------

    >>> from sympy.polys.algebratools import factorial
    >>> factorial(4)
    24
    """
    if not m:
        return 1

    k = m

    while m > 1:
        m -= 1
        k *= m

    return k

class Algebra(object):
    """Represents an abstract domain. """

    dtype = None
    zero  = None
    one   = None

    has_Ring  = False
    has_Field = False

    has_assoc_Ring  = False
    has_assoc_Field = False

    is_ZZ = False
    is_QQ = False

    is_FF = False
    is_CC = False

    is_Poly = False
    is_Frac = False

    is_Exact = True

    is_Numerical = False
    is_Algebraic = False
    is_Composite = False

    has_CharacteristicZero = False

    is_EX = False

    rep = None

    def __init__(self):
        raise NotImplementedError

    def __str__(self):
        return self.rep

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype))

    def __call__(self, *args):
        """Construct an element of `self` domain from `args`. """
        return self.dtype(*args)

    def normal(self, *args):
        """
        Normalize `a` with respect to `self`.

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.polys.algebratools import ZZ_python, ZZ_sympy
        >>> ZZ_python().normal(1)
        1
        >>> type(ZZ_python().normal(1))
        <type 'int'>
        >>> ZZ_sympy().normal(2)
        2
        >>> type(ZZ_sympy().normal(2))
        <class 'sympy.core.numbers.Integer'>
        """
        return self.dtype(*args)

    def convert(K1, a, K0=None):
        """
        Convert an object `a` from `K0` to `K1`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ, ZZ
        >>> QQ.convert(QQ(2, 1), ZZ) == ZZ(2)
        True
        """
        if K0 is not None:
            _convert = getattr(K1, "from_" + K0.__class__.__name__)

            if _convert is not None:
                result = _convert(a, K0)

                if result is not None:
                    return result

            raise CoercionFailed("can't convert %s of type %s to %s" % (a, K0, K1))
        else:
            try:
                if K1.of_type(a):
                    return a

                if type(a) is int:
                    return K1(a)

                if type(a) is long:
                    return K1(a)

                a = sympify(a)

                if isinstance(a, Basic):
                    return K1.from_sympy(a)
            except (TypeError, ValueError):
                pass

            raise CoercionFailed("can't convert %s to type %s" % (a, K1))

    def of_type(self, a):
        """
        Check if `a` is of type `dtype`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ, QQ
        >>> ZZ.of_type(ZZ(2))
        True
        >>> ZZ.of_type(QQ(3, 2))
        False
        """
        return type(a) is type(self.one)

    def __contains__(self, a):
        """Check if `a` belongs to this domain. """
        try:
            self.convert(a)
        except CoercionFailed:
            return False

        return True

    def to_sympy(self, a):
        """
        Convert `a` to a SymPy object.

        Example
        =======
        >>> from sympy import S
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.to_sympy(ZZ(1)) is S.One
        True
        """
        raise NotImplementedError

    def from_sympy(self, a):
        """
        Convert a SymPy object to `dtype`.

        Example
        =======
        >>> from sympy import S
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.from_sympy(S(2))
        2
        """
        raise NotImplementedError

    def from_ZZ_python(K1, a, K0):
        """
        Convert a Python `int` object to `dtype`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ, ZZ_python
        >>> ZZ.from_ZZ_python(2, ZZ_python())
        2
        """
        return None

    def from_QQ_python(K1, a, K0):
        """
        Convert a Python `Fraction` object to `dtype`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> try:
        ...     import fractions # fractions only exists in >= Python 2.6
        ...     from sympy.polys.algebratools import QQ_python
        ...     a = ZZ.from_QQ_python(fractions.Fraction(2, 1), QQ_python())
        ... except ImportError:
        ...     a = ZZ(2)
        >>> a
        2
        """
        return None

    def from_ZZ_sympy(K1, a, K0):
        """
        Convert a SymPy `Integer` object to `dtype`.

        Example
        =======
        >>> from sympy import S
        >>> from sympy.polys.algebratools import ZZ, ZZ_sympy
        >>> ZZ.from_ZZ_sympy(S(2), ZZ_sympy())
        2
        """
        return None

    def from_QQ_sympy(K1, a, K0):
        """
        Convert a SymPy `Rational` object to `dtype`.


        Example
        =======
        >>> from sympy import S
        >>> from sympy.polys.algebratools import QQ, QQ_sympy
        >>> QQ.from_QQ_sympy(S(3)/2, QQ_sympy())
        3/2
        """
        return None

    def from_ZZ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpz` object to `dtype`.


        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import ZZ_gmpy
        ...     a = ZZ.from_ZZ_gmpy(gmpy.mpz(2), ZZ_gmpy)
        ... except ImportError:
        ...     a = ZZ(2)
        >>> a
        2
        """
        return None

    def from_QQ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpq` object to `dtype`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import QQ_gmpy
        ...     a = QQ.from_QQ_gmpy(gmpy.mpq(3, 2), QQ_gmpy())
        ... except ImportError:
        ...     a = QQ(3, 2)
        >>> a
        3/2
        """
        return None

    def from_RR_sympy(K1, a, K0):
        """
        Convert a SymPy `Real` object to `dtype`.

        Example
        =======
        >>> from sympy import S
        >>> from sympy.polys.algebratools import ZZ, RR_sympy
        >>> ZZ.from_RR_sympy(S(1.0), RR_sympy())
        1
        """
        return None

    def from_RR_mpmath(K1, a, K0):
        """
        Convert a mpmath `mpf` object to `dtype`.

        Example
        =======
        >>> from sympy.mpmath import mpf
        >>> from sympy.polys.algebratools import ZZ, RR_mpmath
        >>> ZZ.from_RR_mpmath(1.0, RR_mpmath())
        1
        """
        return None

    def from_AlgebraicField(K1, a, K0):
        """Convert a `ANP` object to `dtype`. """
        return None

    def from_PolynomialRing(K1, a, K0):
        """Convert a `DMP` object to `dtype`. """
        if a.degree() <= 0:
            return K1.convert(a.LC(), K0.dom)

    def from_FractionField(K1, a, K0):
        """Convert a `DMF` object to `dtype`. """
        return None

    def from_ExpressionDomain(K1, a, K0):
        """Convert a `EX` object to `dtype`. """
        return K1.from_sympy(a.ex)

    def unify(K0, K1, gens=None):
        """
        Returns a maximal domain containg `K0` and `K1`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ, QQ
        >>> from sympy.abc import x
        >>> ZZ[x].unify(QQ)
        QQ[x]
        """
        if gens is not None:
            if (K0.is_Composite and (set(K0.gens) & set(gens))) or (K1.is_Composite and (set(K1.gens) & set(gens))):
                raise UnificationFailed("can't unify %s with %s, given %s generators" % (K0, K1, tuple(gens)))
        elif K0 == K1:
            return K0

        if K0.is_EX:
            return K0
        if K1.is_EX:
            return K1

        if not K0.is_Exact:
            return K0
        if not K1.is_Exact:
            return K1

        if K0.is_Composite:
            if K1.is_Composite:
                if K0.gens == K1.gens:
                    if K0.has_Field and K1.has_Field:
                        if K0.dom.has_Field:
                            return K0
                        else:
                            return K1
                    elif K0.has_Field:
                        if K0.dom == K1.dom:
                            return K0
                    elif K1.has_Field:
                        if K0.dom == K1.dom:
                            return K1
                    else:
                        if K0.dom.has_Field:
                            return K0
                        else:
                            return K1
                else:
                    gens = set(K0.gens + K1.gens)

                    try:
                        gens = sorted(gens)
                    except TypeError:
                        gens = list(gens)

                    if K0.has_Field and K1.has_Field:
                        if K0.dom.has_Field:
                            return K0.__class__(K0.dom, *gens)
                        else:
                            return K1.__class__(K1.dom, *gens)
                    elif K0.has_Field:
                        if K0.dom == K1.dom:
                            return K0.__class__(K0.dom, *gens)
                    elif K1.has_Field:
                        if K0.dom == K1.dom:
                            return K1.__class__(K1.dom, *gens)
                    else:
                        if K0.dom.has_Field:
                            return K0.__class__(K0.dom, *gens)
                        else:
                            return K1.__class__(K1.dom, *gens)
            elif K1.is_Algebraic:
                raise UnificationFailed("can't unify %s with %s" % (K0, K1))
            else:
                if K0.has_Field:
                    if K0.dom == K1:
                        return K0
                else:
                    if K0.dom.has_Field:
                        return K0
                    else:
                        return K0.__class__(K1, *K0.gens)
        elif K0.is_Algebraic:
            if K1.is_Composite:
                raise UnificationFailed("can't unify %s with %s" % (K0, K1))
            elif K1.is_Algebraic:
                raise NotImplementedError("unification of different algebraic extensions")
            elif K1.is_ZZ or K1.is_QQ:
                return K0
            else:
                raise UnificationFailed("can't unify %s with %s" % (K0, K1))
        else:
            if K1.is_Composite:
                if K1.has_Field:
                    if K0 == K1.dom:
                        return K1
                else:
                    if K1.dom.has_Field:
                        return K1
                    else:
                        return K1.__class__(K0, *K1.gens)
            elif K1.is_Algebraic:
                if K0.is_ZZ or K0.is_QQ:
                    return K1
                else:
                    raise UnificationFailed("can't unify %s with %s" % (K0, K1))
            else:
                if K0.has_Field:
                    return K0
                else:
                    return K1

        return ExpressionDomain()

    def __eq__(self, other):
        """Returns `True` if two algebras are equivalent. """
        return self.dtype == other.dtype

    def __ne__(self, other):
        """Returns `False` if two algebras are equivalent. """
        return self.dtype != other.dtype

    def get_ring(self):
        """
        Returns a ring associated with `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.get_ring()
        ZZ
        """
        raise DomainError('there is no ring associated with %s' % self)

    def get_field(self):
        """
        Returns a field associated with `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.get_field()
        QQ
        """
        raise DomainError('there is no field associated with %s' % self)

    def get_exact(self):
        """
        Returns an exact domain associated with `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import RR
        >>> RR.get_exact()
        QQ
        """
        return self

    def float_domain(self):
        """
        Returns the floating point domain associated with `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.float_domain()
        FF
        """
        return FF

    def complex_domain(self):
        """
        Returns the complex domain associated with `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.complex_domain()
        CC
        """
        return CC

    def __getitem__(self, gens):
        """The mathematical way do make a polynomial ring. """
        if hasattr(gens, '__iter__'):
            return self.poly_ring(*gens)
        else:
            return self.poly_ring(gens)

    def poly_ring(self, *gens):
        """
        Returns a polynomial ring, i.e. `K[X]`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> from sympy.abc import x
        >>> ZZ.poly_ring(x)
        ZZ[x]
        >>> ZZ[x]
        ZZ[x]
        """
        return PolynomialRing(self, *gens)

    def frac_field(self, *gens):
        """
        Returns a fraction field, i.e. `K(X)`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> from sympy.abc import x
        >>> ZZ.frac_field(x)
        ZZ(x)
        """
        return FractionField(self, *gens)

    def algebraic_field(self, *extension):
        """
        Returns an algebraic field, i.e. `K(alpha, ...)`.

        Example
        =======
        >>> from sympy import I
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.algebraic_field(I)
        QQ<I>
        """
        raise DomainError("can't create algebraic field over %s" % self)

    def is_zero(self, a):
        """
        Returns True if `a` is zero.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.is_zero(ZZ(0))
        True
        >>> ZZ.is_zero(ZZ(1))
        False
        """
        return not a

    def is_one(self, a):
        """
        Returns True if `a` is one.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.is_one(ZZ(1))
        True
        >>> ZZ.is_one(ZZ(0))
        False
        """
        return a == self.one

    def is_positive(self, a):
        """
        Returns True if `a` is positive.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.is_positive(ZZ(2))
        True
        >>> ZZ.is_positive(ZZ(0))
        False
        >>> ZZ.is_positive(ZZ(-2))
        False
        """
        return a > 0

    def is_negative(self, a):
        """
        Returns True if `a` is negative.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.is_negative(ZZ(2))
        False
        >>> ZZ.is_negative(ZZ(0))
        False
        >>> ZZ.is_negative(ZZ(-2))
        True
        """
        return a < 0

    def is_nonpositive(self, a):
        """
        Returns True if `a` is non-positive.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.is_nonpositive(ZZ(2))
        False
        >>> ZZ.is_nonpositive(ZZ(0))
        True
        >>> ZZ.is_nonpositive(ZZ(-2))
        True
        """
        return a <= 0

    def is_nonnegative(self, a):
        """
        Returns True if `a` is non-negative.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.is_nonnegative(ZZ(2))
        True
        >>> ZZ.is_nonnegative(ZZ(0))
        True
        >>> ZZ.is_nonnegative(ZZ(-2))
        False
        """
        return a >= 0

    def abs(self, a):
        """
        Absolute value of `a`, implies `__abs__`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.abs(ZZ(-2))
        2
        """
        return abs(a)

    def neg(self, a):
        """
        Returns `a` negated, implies `__neg__`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.neg(ZZ(2))
        -2
        """
        return -a

    def pos(self, a):
        """
        Returns `a` positive, implies `__pos__`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.pos(ZZ(2))
        2
        >>> ZZ.pos(ZZ(-2))
        -2
        """
        return +a

    def add(self, a, b):
        """
        Sum of `a` and `b`, implies `__add__`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.add(ZZ(1), ZZ(2))
        3
        """
        return a + b

    def sub(self, a, b):
        """
        Difference of `a` and `b`, implies `__sub__`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.sub(ZZ(1), ZZ(2))
        -1
        """
        return a - b

    def mul(self, a, b):
        """
        Product of `a` and `b`, implies `__mul__`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.mul(ZZ(2), ZZ(3))
        6
        """
        return a * b

    def pow(self, a, b):
        """
        Raise `a` to power `b`, implies `__pow__`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.pow(ZZ(2), ZZ(3))
        8
        """
        return a ** b

    def exquo(self, a, b):
        """
        Exact quotient of `a` and `b`, implies something.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.exquo(ZZ(3), ZZ(2))
        1
        >>> ZZ.exquo(ZZ(4), ZZ(2))
        2
        """
        raise NotImplementedError

    def quo(self, a, b):
        """
        Quotient of `a` and `b`, implies something.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.quo(ZZ(5), ZZ(2))
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: 2 does not divide 5 in ZZ
        >>> ZZ.quo(ZZ(4), ZZ(2))
        2
        """
        raise NotImplementedError

    def rem(self, a, b):
        """
        Remainder of `a` and `b`, implies `__mod__`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.rem(ZZ(5), ZZ(2))
        1
        """
        raise NotImplementedError

    def div(self, a, b):
        """
        Division of `a` and `b`, implies something.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.div(ZZ(5), ZZ(2))
        (2, 1)
        """
        raise NotImplementedError

    def invert(self, a, b):
        """
        Returns the inversion of `a mod b`, implies something.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.invert(ZZ(3), ZZ(5))
        2
        """
        raise NotImplementedError

    def numer(self, a):
        """
        Returns the numerator of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.numer(QQ(3, 2))
        3
        """
        raise NotImplementedError

    def denom(self, a):
        """
        Returns the denominator of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.denom(QQ(3, 2))
        2
        """
        raise NotImplementedError

    def gcdex(self, a, b):
        """
        Extended GCD of `a` and `b`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> tuple(ZZ.gcdex(ZZ(12), ZZ(8)))
        (1, -1, 4)
        """
        raise NotImplementedError

    def gcd(self, a, b):
        """
        Returns the GCD of `a` and `b`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.gcd(ZZ(12), ZZ(8))
        4
        """
        raise NotImplementedError

    def lcm(self, a, b):
        """
        Returns the LCM of `a` and `b`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.lcm(ZZ(12), ZZ(8))
        24
        """
        raise NotImplementedError

    def log(self, a, b):
        """
        Returns the b-base logarithm of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.log(ZZ(8), 2)
        3
        """
        raise NotImplementedError

    def sqrt(self, a):
        """
        Returns the square root of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.sqrt(ZZ(9))
        3
        """
        raise NotImplementedError

    def factorial(self, a):
        """
        Returns the factorial of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.factorial(4)
        24
        """
        return self.dtype(factorial(a))

    def evalf(self, a, **args):
        """
        Returns a numerical approximation of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.evalf(QQ(3, 2))
        1.50000000000000
        """
        return self.to_sympy(a).evalf(**args)

    def real(self, a):
        """
        Returns the real part of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import CC
        >>> CC.real(CC(1+2j))
        1.0
        """
        return a

    def imag(self, a):
        """
        Returns the imaginary part of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import CC
        >>> CC.imag(CC(1+2j))
        2.0
        """
        return self.zero

class Ring(Algebra):
    """Represents a ring domain. """

    has_Ring = True

    def get_ring(self):
        """
        Returns a ring associated with `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.get_ring()
        ZZ
        """
        return self

    def exquo(self, a, b):
        """
        Exact quotient of `a` and `b`, implies `__floordiv__`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.exquo(ZZ(3), ZZ(2))
        1
        """
        return a // b

    def quo(self, a, b):
        """
        Quotient of `a` and `b`, implies `__floordiv__`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.quo(ZZ(3), ZZ(2))
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: 2 does not divide 3 in ZZ
        >>> ZZ.quo(ZZ(4), ZZ(2))
        2
        """
        if a % b:
            raise ExactQuotientFailed('%s does not divide %s in %s' % (b, a, self))
        else:
            return a // b

    def rem(self, a, b):
        """
        Remainder of `a` and `b`, implies `__mod__`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.rem(ZZ(3), ZZ(2))
        1
        """
        return a % b

    def div(self, a, b):
        """
        Division of `a` and `b`, implies `__divmod__`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.div(ZZ(3), ZZ(2))
        (1, 1)
        """
        return divmod(a, b)

    def invert(self, a, b):
        """
        Returns the inversion of `a mod b`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.invert(ZZ(3), ZZ(5))
        2
        """
        s, t, h = self.gcdex(a, b)

        if self.is_one(h):
            return s % b
        else:
            raise NotInvertible("zero divisor")

    def numer(self, a):
        """
        Returns the numerator of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.numer(ZZ(2))
        2
        """
        return a

    def denom(self, a):
        """
        Returns the denominator of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.denom(ZZ(2))
        1
        """
        return self.one

class Field(Ring):
    """Represents a field domain. """

    has_Field = True

    def get_field(self):
        """
        Returns a field associated with `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.get_field()
        QQ
        """
        return self

    def exquo(self, a, b):
        """
        Exact quotient of `a` and `b`, implies `__div__`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.exquo(QQ(3), QQ(2))
        3/2
        """
        return a / b

    def quo(self, a, b):
        """
        Quotient of `a` and `b`, implies `__div__`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.quo(QQ(3), QQ(2))
        3/2
        """
        return a / b

    def rem(self, a, b):
        """
        Remainder of `a` and `b`, implies nothing.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.rem(QQ(3), QQ(2)) == QQ(0)
        True
        """
        return self.zero

    def div(self, a, b):
        """
        Division of `a` and `b`, implies `__div__`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.div(QQ(3), QQ(2)) == (QQ(3, 2), QQ(0))
        True
        """
        return a / b, self.zero

    def gcd(self, a, b):
        """
        Returns the GCD of `a` and `b`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.gcd(QQ(3, 2), QQ(5, 6)) == QQ(1)
        True
        """
        return self.one

    def lcm(self, a, b):
        """
        Returns the LCM of `a` and `b`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.lcm(QQ(3, 2), QQ(5, 6))
        5/4
        """
        return a*b

class IntegerRing(Ring):
    """General class for integer rings. """

    is_ZZ = True
    rep   = 'ZZ'

    is_Numerical = True

    has_assoc_Ring         = True
    has_assoc_Field        = True

    has_CharacteristicZero = True

    def get_field(self):
        """
        Returns a field associated with `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.get_field()
        QQ
        """
        return QQ

    def from_AlgebraicField(K1, a, K0):
        """Convert a `ANP` object to `dtype`. """
        if a.is_ground:
            return K1.convert(a.LC(), K0.dom)

    def log(self, a, b):
        """
        Returns the b-base logarithm of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.log(8, 2)
        3
        """
        return self.dtype(math.log(a, b))

class RationalField(Field):
    """General class for rational fields. """

    is_QQ = True
    rep   = 'QQ'

    is_Numerical = True

    has_assoc_Ring         = True
    has_assoc_Field        = True

    has_CharacteristicZero = True

    def get_ring(self):
        """
        Returns a ring associated with `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.get_ring()
        ZZ
        """
        return ZZ

    def algebraic_field(self, *extension):
        """
        Returns an algebraic field, i.e. `QQ(alpha, ...)`.

        Example
        =======
        >>> from sympy import sqrt
        >>> from sympy.polys.algebratools import QQ
        >>> QQ.algebraic_field(sqrt(2) + sqrt(3))
        QQ<2**(1/2) + 3**(1/2)>
        """
        return AlgebraicField(self, *extension)

    def from_AlgebraicField(K1, a, K0):
        """Convert a `ANP` object to `dtype`. """
        if a.is_ground:
            return K1.convert(a.LC(), K0.dom)

class ZZ_python(IntegerRing):
    pass
class QQ_python(RationalField):
    pass
class ZZ_sympy(IntegerRing):
    pass
class QQ_sympy(RationalField):
    pass
class ZZ_gmpy(IntegerRing):
    pass
class QQ_gmpy(RationalField):
    pass

class PolynomialRing(Ring):
    pass
class FractionField(Field):
    pass

class ExpressionDomain(Field):
    pass

HAS_FRACTION = True

try:
    import fractions
except ImportError:
    HAS_FRACTION = False

HAS_GMPY = True

try:
    import gmpy
except ImportError:
    HAS_GMPY = False

from sympy.core.numbers import (
    igcdex as python_gcdex,
    igcd as python_gcd,
    ilcm as python_lcm,
)

from sympy.mpmath.libmp.libmpf import (
    isqrt as python_sqrt,
)

from __builtin__ import (
    int as python_int,
)

if HAS_FRACTION:
    from fractions import (
        Fraction as python_rat,
    )

class ZZ_python(IntegerRing):
    """Integer ring based on Python int class. """

    dtype = python_int
    zero  = dtype(0)
    one   = dtype(1)

    def __init__(self):
        pass

    def to_sympy(self, a):
        """
        Convert `a` to a SymPy object.

        Example
        =======
        >>> from sympy import S
        >>> from sympy.polys.algebratools import ZZ_python
        >>> ZZ_python().to_sympy(1)
        1
        >>> type(ZZ_python().to_sympy(1))
        <class 'sympy.core.numbers.One'>
        """
        return sympy_int(a)

    def from_sympy(self, a):
        """
        Convert SymPy's Integer to `dtype` (`int`).

        Example
        =======
        >>> from sympy import S
        >>> from sympy.polys.algebratools import ZZ_python
        >>> ZZ_python().from_sympy(S.One)
        1
        >>> type(ZZ_python().from_sympy(S.One))
        <type 'int'>
        """
        if a.is_Integer:
            return python_int(a.p)
        elif a.is_Real and int(a) == a:
            return sympy_int(int(a))
        else:
            raise CoercionFailed("expected `Integer` object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """
        Convert a Python `int` object to `dtype` (`int`).

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ_python
        >>> ZZ_python().from_ZZ_python(int(1), ZZ_python())
        1
        >>> type(ZZ_python().from_ZZ_python(int(1), ZZ_python()))
        <type 'int'>
        """
        return a

    def from_QQ_python(K1, a, K0):
        """
        Convert a Python `Fraction` object to `dtype` (`int`).

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ_python
        >>> try:
        ...     import fractions # fractions only exists in >= Python 2.6
        ...     from sympy.polys.algebratools import QQ_python
        ...     a = ZZ_python().from_QQ_python(fractions.Fraction(3, 1),
        ...     QQ_python())
        ... except ImportError:
        ...     a = 3
        >>> a
        3
        >>> type(a)
        <type 'int'>
        """
        # XXX: is there a better way to structure this doctest?
        if a.denominator == 1:
            return a.numerator

    def from_ZZ_sympy(K1, a, K0):
        """
        Convert a SymPy `Integer` object to `dtype` (`int`).

        Example
        =======
        >>> from sympy import S
        >>> from sympy.polys.algebratools import ZZ_python, ZZ_sympy
        >>> ZZ_python().from_ZZ_sympy(S.One, ZZ_sympy())
        1
        >>> type(ZZ_python().from_ZZ_sympy(S.One, ZZ_sympy()))
        <type 'int'>
        """
        return a.p

    def from_QQ_sympy(K1, a, K0):
        """
        Convert a SymPy `Rational` object to `dtype` (`int`).

        Example
        =======
        >>> from sympy import Rational
        >>> from sympy.polys.algebratools import ZZ_python, QQ_sympy
        >>> ZZ_python().from_QQ_sympy(Rational(2, 1), QQ_sympy())
        2
        >>> type(ZZ_python().from_QQ_sympy(Rational(2, 1), QQ_sympy()))
        <type 'int'>
        """
        if a.q == 1:
            return a.p

    def from_ZZ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpz` object to `dtype` (`int).

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ_python
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import ZZ_gmpy
        ...     a = ZZ_python().from_ZZ_gmpy(gmpy.mpz(1), ZZ_gmpy())
        ... except ImportError:
        ...     a = 1
        >>> a
        1
        >>> type(a)
        <type 'int'>
        """
        return python_int(a)

    def from_QQ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpq` object to `dtype` (`int).

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ_python
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import QQ_gmpy
        ...     a = ZZ_python().from_QQ_gmpy(gmpy.mpq(2, 1), QQ_gmpy())
        ... except ImportError:
        ...     a = 2
        >>> a
        2
        >>> type(a)
        <type 'int'>
        """
        if a.denom() == 1:
            return python_int(a.numer())

    def from_RR_sympy(K1, a, K0):
        """
        Convert a SymPy `Real` object to `dtype` (`int`).

        Example
        =======
        >>> from sympy import Real
        >>> from sympy.polys.algebratools import ZZ_python, RR_sympy
        >>> ZZ_python().from_RR_sympy(Real(1.0), RR_sympy())
        1
        >>> type(ZZ_python().from_RR_sympy(Real(1.0), RR_sympy()))
        <type 'int'>
        """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return python_int(p)

    def from_RR_mpmath(K1, a, K0):
        """
        Convert a mpmath `mpf` object to `dtype` (`int`).

        Example
        =======
        >>> from sympy.mpmath import mpf
        >>> from sympy.polys.algebratools import ZZ_python, RR_mpmath
        >>> ZZ_python().from_RR_mpmath(mpf(1), RR_mpmath())
        1
        >>> type(ZZ_python().from_RR_mpmath(mpf(1), RR_mpmath()))
        <type 'int'>
        """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return python_int(p)

    def gcdex(self, a, b):
        """
        Extended GCD of `a` and `b`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ_python
        >>> ZZ_python().gcdex(12, 8)
        (1, -1, 4)
        >>> map(type, ZZ_python().gcdex(12, 8))
        [<type 'int'>, <type 'int'>, <type 'int'>]
        """
        return python_gcdex(a, b)

    def gcd(self, a, b):
        """
        Returns the GCD of `a` and `b`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ_python
        >>> ZZ_python().gcd(12, 8)
        4
        >>> type(ZZ_python().gcd(12, 8))
        <type 'int'>
        """
        return python_gcd(a, b)

    def lcm(self, a, b):
        """
        Returns the LCM of `a` and `b`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ_python
        >>> ZZ_python().lcm(12, 8)
        24
        >>> type(ZZ_python().lcm(12, 8))
        <type 'int'>
        """
        return python_lcm(a, b)

    def sqrt(self, a):
        """
        Returns the square root of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ_python
        >>> ZZ_python().sqrt(9)
        3
        >>> ZZ_python().sqrt(12)
        3
        >>> type(ZZ_python().sqrt(9))
        <type 'int'>
        """
        return python_int(python_sqrt(a))

if HAS_FRACTION:
    class QQ_python(RationalField):
        """Rational field based on Python Fraction class. """

        dtype = python_rat
        zero  = dtype(0)
        one   = dtype(1)

        def __init__(self):
            pass

        def to_sympy(self, a):
            """
            Convert `a` to a SymPy object.

            Example
            =======
            >>> from sympy import Rational
            >>> import fractions
            >>> from sympy.polys.algebratools import QQ_python
            >>> QQ_python().to_sympy(fractions.Fraction(2, 3))
            2/3
            >>> type(QQ_python().to_sympy(fractions.Fraction(2, 3)))
            <class 'sympy.core.numbers.Rational'>
            """
            return sympy_rat(a.numerator, a.denominator)

        def from_sympy(self, a):
            """
            Convert SymPy's Rational to `dtype` (`fractions.Fraction`).

            Example
            =======
            >>> from sympy import Rational
            >>> from sympy.polys.algebratools import QQ_python
            >>> QQ_python().from_sympy(Rational(2, 3))
            2/3
            >>> type(QQ_python().from_sympy(Rational(2, 3)))
            <class 'fractions.Fraction'>
            """
            if a.is_Rational and a.q != 0:
                return python_rat(a.p, a.q)
            elif a.is_Real:
                return python_rat(*RR.as_integer_ratio(a))
            else:
                raise CoercionFailed("expected `Rational` object, got %s" % a)

        def from_ZZ_python(K1, a, K0):
            """
            Convert a Python `int` object to `dtype` (`fractions.Fraction`).

            Example
            =======
            >>> from sympy.polys.algebratools import QQ_python, ZZ_python
            >>> QQ_python().from_ZZ_python(1, ZZ_python())
            1/1
            >>> type(QQ_python().from_ZZ_python(1, ZZ_python()))
            <class 'fractions.Fraction'>
            """
            return python_rat(a)

        def from_QQ_python(K1, a, K0):
            """
            Convert a Python `Fraction` object to `dtype` (`fractions.Fraction`).

            Example
            =======
            >>> import fractions
            >>> from sympy.polys.algebratools import QQ_python
            >>> QQ_python().from_QQ_python(fractions.Fraction(2, 3), QQ_python())
            2/3
            >>> type(QQ_python().from_QQ_python(fractions.Fraction(2, 3),
            ... QQ_python()))
            <class 'fractions.Fraction'>
            """
            return a

        def from_ZZ_sympy(K1, a, K0):
            """
            Convert a SymPy `Integer` object to `dtype` (`fractions.Fraction`).

            Example
            =======
            >>> from sympy import S
            >>> from sympy.polys.algebratools import QQ_python, ZZ_sympy
            >>> QQ_python().from_ZZ_sympy(S.One, ZZ_sympy())
            1/1
            >>> type(QQ_python().from_ZZ_sympy(S.One, ZZ_sympy()))
            <class 'fractions.Fraction'>
            """
            return python_rat(a.p)

        def from_QQ_sympy(K1, a, K0):
            """
            Convert a SymPy `Rational` object to `dtype` (`fractions.Fraction`).

            Example
            =======
            >>> from sympy import Rational
            >>> from sympy.polys.algebratools import QQ_python, QQ_sympy
            >>> QQ_python().from_QQ_sympy(Rational(2, 3), QQ_sympy())
            2/3
            >>> type(QQ_python().from_QQ_sympy(Rational(2, 3), QQ_sympy()))
            <class 'fractions.Fraction'>
            """
            return python_rat(a.p, a.q)

        def from_ZZ_gmpy(K1, a, K0):
            """
            Convert a GMPY `mpz` object to `dtype` (`fractions.Fraction`).

            Example
            =======
            >>> import fractions
            >>> try:
            ...     import gmpy
            ...     from sympy.polys.algebratools import QQ_python, ZZ_gmpy
            ...     a = QQ_python().from_ZZ_gmpy(gmpy.mpz(1), ZZ_gmpy())
            ... except ImportError:
            ...     a = fractions.Fraction(1, 1)
            >>> a
            1/1
            >>> type(a)
            <class 'fractions.Fraction'>
            """
            return python_rat(python_int(a))

        def from_QQ_gmpy(K1, a, K0):
            """
            Convert a GMPY `mpq` object to `dtype`(`fractions.Fraction`).

            Example
            =======
            >>> import fractions
            >>> try:
            ...     import gmpy
            ...     from sympy.polys.algebratools import QQ_python, QQ_gmpy
            ...     a = QQ_python().from_QQ_gmpy(gmpy.mpq(2, 3), QQ_gmpy())
            ... except ImportError:
            ...     a = fractions.Fraction(2, 3)
            >>> a
            2/3
            >>> type(a)
            <class 'fractions.Fraction'>
            """
            return python_rat(python_int(a.numer()),
                              python_int(a.denom()))

        def from_RR_sympy(K1, a, K0):
            """
            Convert a SymPy `Real` object to `dtype` (`fractions.Fraction`).

            Example
            =======
            >>> from sympy import Real
            >>> from sympy.polys.algebratools import QQ_python, RR_sympy
            >>> QQ_python().from_RR_sympy(Real(1.5), RR_sympy())
            3/2
            >>> type(QQ_python().from_RR_sympy(Real(1.5), RR_sympy()))
            <class 'fractions.Fraction'>
            """
            return python_rat(*K0.as_integer_ratio(a))

        def from_RR_mpmath(K1, a, K0):
            """
            Convert a mpmath `mpf` object to `dtype` (`fractions.Fraction`).

            Example
            =======
            >>> from sympy.mpmath import mpf
            >>> from sympy.polys.algebratools import QQ_python, RR_mpmath
            >>> QQ_python().from_RR_mpmath(mpf(1.5), RR_mpmath())
            3/2
            >>> type(QQ_python().from_RR_mpmath(mpf(1.5), RR_mpmath()))
            <class 'fractions.Fraction'>
            """
            return python_rat(*K0.as_integer_ratio(a))

        def numer(self, a):
            """
            Returns the numerator of `a`.

            Example
            =======
            >>> import fractions
            >>> from sympy.polys.algebratools import QQ_python
            >>> QQ_python().numer(fractions.Fraction(2, 3))
            2
            >>> type(QQ_python().numer(fractions.Fraction(2, 3)))
            <type 'int'>
            """
            return a.numerator

        def denom(self, a):
            """
            Returns the denominator of `a`.

            Example
            =======
            >>> import fractions
            >>> from sympy.polys.algebratools import QQ_python
            >>> QQ_python().denom(fractions.Fraction(2, 3))
            3
            >>> type(QQ_python().denom(fractions.Fraction(2, 3)))
            <type 'int'>
            """
            return a.denominator

from sympy import (
    Integer as sympy_int,
    Rational as sympy_rat,
)

from sympy.core.numbers import (
    igcdex as sympy_gcdex,
    igcd as sympy_gcd,
    ilcm as sympy_lcm,
)

from sympy.mpmath.libmp.libmpf import (
    isqrt as sympy_sqrt,
)

class ZZ_sympy(IntegerRing):
    """Integer ring based on SymPy Integer class. """

    dtype = sympy_int
    zero  = dtype(0)
    one   = dtype(1)

    def __init__(self):
        pass

    def of_type(self, a):
        """
        Check if `a` is of type `dtype` (`sympy`).

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.polys.algebratools import ZZ_sympy, QQ
        >>> ZZ_sympy().of_type(Integer(2))
        True
        >>> ZZ_sympy().of_type(QQ(3, 2))
        False
        """
        return type(a) in [type(self.one), type(self.zero),
                           type(sympy_int(-1)), type(sympy_int(2))]


    def to_sympy(self, a):
        """
        Convert `a` to a SymPy object.

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.polys.algebratools import ZZ_sympy
        >>> ZZ_sympy().to_sympy(Integer(2))
        2
        >>> type(ZZ_sympy().to_sympy(Integer(2)))
        <class 'sympy.core.numbers.Integer'>
        """
        return a

    def from_sympy(self, a):
        """
        Convert SymPy's Integer to `dtype` (`Integer`).

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.polys.algebratools import ZZ_sympy
        >>> ZZ_sympy().from_sympy(Integer(2))
        2
        >>> type(ZZ_sympy().from_sympy(Integer(2)))
        <class 'sympy.core.numbers.Integer'>
        """
        if a.is_Integer:
            return a
        elif a.is_Real and int(a) == a:
            return sympy_int(int(a))
        else:
            raise CoercionFailed("expected Integer object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """
        Convert a Python `int` object to `dtype` (`Integer`).

        Example
        =======
        >>> from sympy.polys.algebratools import ZZ_sympy, ZZ_python
        >>> ZZ_sympy().from_ZZ_python(2, ZZ_python())
        2
        >>> type(ZZ_sympy().from_ZZ_python(2, ZZ_python()))
        <class 'sympy.core.numbers.Integer'>
        """
        return sympy_int(a)

    def from_QQ_python(K1, a, K0):
        """
        Convert a Python `Fraction` object to `dtype` (`Integer`).

        Example
        =======
        >>> from sympy import Integer
        >>> try:
        ...     import fractions # fractions only exists in >= Python 2.6
        ...     from sympy.polys.algebratools import ZZ_sympy, QQ_python
        ...     a = ZZ_sympy().from_QQ_python(fractions.Fraction(2, 1),
        ...     QQ_python())
        ... except ImportError:
        ...     a = Integer(2)
        >>> a
        2
        >>> type(a)
        <class 'sympy.core.numbers.Integer'>
        """
        if a.denominator == 1:
            return sympy_int(a.numerator)

    def from_ZZ_sympy(K1, a, K0):
        """
        Convert a SymPy `Integer` object to `dtype` (`Integer`).

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.polys.algebratools import ZZ_sympy
        >>> ZZ_sympy().from_ZZ_sympy(Integer(2), ZZ_sympy())
        2
        >>> type(ZZ_sympy().from_ZZ_sympy(Integer(2), ZZ_sympy()))
        <class 'sympy.core.numbers.Integer'>
        """
        return a

    def from_QQ_sympy(K1, a, K0):
        """
        Convert a SymPy `Rational` object to `dtype` (`Integer`).

        Example
        =======
        >>> from sympy import Rational
        >>> from sympy.polys.algebratools import ZZ_sympy, QQ_sympy
        >>> ZZ_sympy().from_QQ_sympy(Rational(2, 1), QQ_sympy())
        2
        >>> type(ZZ_sympy().from_QQ_sympy(Rational(2, 1), QQ_sympy()))
        <class 'sympy.core.numbers.Integer'>
        """
        if a.q == 1:
            return sympy_int(a.p)

    def from_ZZ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpz` object to `dtype` (`Integer`).

        Example
        =======
        >>> from sympy import Integer
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import ZZ_sympy, ZZ_gmpy
        ...     a = ZZ_sympy().from_ZZ_gmpy(gmpy.mpz(2), ZZ_gmpy())
        ... except ImportError:
        ...     a = Integer(2)
        >>> a
        2
        >>> type(a)
        <class 'sympy.core.numbers.Integer'>
        """
        return sympy_int(python_int(a))

    def from_QQ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpq` object to `dtype` (`Integer`).

        Example
        =======
        >>> from sympy import Integer
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import ZZ_sympy, QQ_gmpy
        ...     a = ZZ_sympy().from_QQ_gmpy(gmpy.mpq(2, 1), QQ_gmpy())
        ... except ImportError:
        ...     a = Integer(2)
        >>> a
        2
        >>> type(a)
        <class 'sympy.core.numbers.Integer'>
        """
        if a.denom() == 1:
            return sympy_int(python_int(a.numer()))

    def from_RR_sympy(K1, a, K0):
        """
        Convert a SymPy `Real` object to `dtype` (`Integer`).

        Example
        =======
        >>> from sympy import Real
        >>> from sympy.polys.algebratools import ZZ_sympy, RR_sympy
        >>> ZZ_sympy().from_RR_sympy(Real(2.0), RR_sympy())
        2
        >>> type(ZZ_sympy().from_RR_sympy(Real(2.0), RR_sympy()))
        <class 'sympy.core.numbers.Integer'>
        """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return sympy_int(p)

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype` (`Integer`).

        Example
        =======
        >>> from sympy.mpmath import mpf
        >>> from sympy.polys.algebratools import ZZ_sympy, RR_mpmath
        >>> ZZ_sympy().from_RR_mpmath(mpf(2.0), RR_mpmath())
        2
        >>> type(ZZ_sympy().from_RR_mpmath(mpf(2.0), RR_mpmath()))
        <class 'sympy.core.numbers.Integer'>
        """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return sympy_int(p)

    def gcdex(self, a, b):
        """
        Extended GCD of `a` and `b`.

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.polys.algebratools import ZZ_sympy
        >>> ZZ_sympy().gcdex(Integer(12), Integer(8))
        [1, -1, 4]
        >>> map(type, ZZ_sympy().gcdex(Integer(12), Integer(8)))
        [<class 'sympy.core.numbers.One'>,
         <class 'sympy.core.numbers.NegativeOne'>,
         <class 'sympy.core.numbers.Integer'>]
        """
        return map(sympy_int, sympy_gcdex(int(a), int(b)))

    def gcd(self, a, b):
        """
        Returns the GCD of `a` and `b`.

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.polys.algebratools import ZZ_sympy
        >>> ZZ_sympy().gcd(Integer(12), Integer(8))
        4
        >>> type(ZZ_sympy().gcd(Integer(12), Integer(8)))
        <class 'sympy.core.numbers.Integer'>
        """
        return sympy_int(sympy_gcd(int(a), int(b)))

    def lcm(self, a, b):
        """
        Returns the LCM of `a` and `b`.

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.polys.algebratools import ZZ_sympy
        >>> ZZ_sympy().lcm(Integer(12), Integer(8))
        24
        >>> type(ZZ_sympy().lcm(Integer(12), Integer(8)))
        <class 'sympy.core.numbers.Integer'>
        """
        return sympy_int(sympy_lcm(int(a), int(b)))

    def sqrt(self, a):
        """
        Returns the square root of `a`.

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.polys.algebratools import ZZ_sympy
        >>> ZZ_sympy().sqrt(Integer(9))
        3
        >>> ZZ_sympy().sqrt(Integer(12))
        3
        >>> type(ZZ_sympy().sqrt(Integer(9)))
        <class 'sympy.core.numbers.Integer'>
        """
        return sympy_int(int(sympy_sqrt(int(a))))

class QQ_sympy(RationalField):
    """Rational field based on SymPy Rational class. """

    dtype = sympy_rat
    zero  = dtype(0)
    one   = dtype(1)

    def __init__(self):
        pass

    def of_type(self, a):
        """
        Check if `a` is of type `dtype` (`Rational`).

        Example
        =======
        >>> from sympy import Rational, Real
        >>> from sympy.polys.algebratools import QQ_sympy
        >>> QQ_sympy().of_type(Rational(3, 2))
        True
        >>> QQ_sympy().of_type(2)
        False
        """
        return type(a) in [type(self.one), type(self.zero), type(sympy_rat(-1)),
                           type(sympy_rat(2)), type(sympy_rat(1, 2)),
                           type(sympy_rat(3, 2))]


    def to_sympy(self, a):
        """
        Convert `a` to a SymPy object.

        Example
        =======
        >>> from sympy import Rational
        >>> from sympy.polys.algebratools import QQ_sympy
        >>> QQ_sympy().to_sympy(Rational(3, 2))
        3/2
        >>> type(QQ_sympy().to_sympy(Rational(3, 2)))
        <class 'sympy.core.numbers.Rational'>
        """
        return a

    def from_sympy(self, a):
        """
        Convert SymPy's Rational to `dtype` (`Rational`).

        Example
        =======
        >>> from sympy import Rational
        >>> from sympy.polys.algebratools import QQ_sympy
        >>> QQ_sympy().from_sympy(Rational(3, 2))
        3/2
        >>> type(QQ_sympy().from_sympy(Rational(3, 2)))
        <class 'sympy.core.numbers.Rational'>
        """
        if a.is_Rational and a.q != 0:
            return a
        elif a.is_Real:
            return sympy_rat(*RR.as_integer_ratio(a))
        else:
            raise CoercionFailed("expected `Rational` object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """
        Convert a Python `int` object to `dtype` (`Rational`).

        Example
        =======
        >>> from sympy.polys.algebratools import QQ_sympy, ZZ_python
        >>> QQ_sympy().from_ZZ_python(2, ZZ_python())
        2
        >>> type(QQ_sympy().from_ZZ_python(2, ZZ_python()))
        <class 'sympy.core.numbers.Integer'>
        """
        return sympy_rat(a)

    def from_QQ_python(K1, a, K0):
        """
        Convert a Python `Fraction` object to `dtype` (`Rational`).

        Example
        =======
        >>> from sympy import Rational
        >>> try:
        ...     import fractions # fractions only exists in >= Python 2.6
        ...     from sympy.polys.algebratools import QQ_sympy, QQ_python
        ...     a = QQ_sympy().from_QQ_python(fractions.Fraction(3, 2),
        ...     QQ_python())
        ... except ImportError:
        ...     a = Rational(3, 2)
        >>> a
        3/2
        >>> type(a)
        <class 'sympy.core.numbers.Rational'>
        """
        return sympy_rat(a.numerator, a.denominator)

    def from_ZZ_sympy(K1, a, K0):
        """
        Convert a SymPy `Integer` object to `dtype` (`Rational`).

        Example
        =======
        >>> from sympy import Rational, Integer
        >>> from sympy.polys.algebratools import QQ_sympy, ZZ_sympy
        >>> QQ_sympy().from_ZZ_sympy(Integer(2), ZZ_sympy())
        2
        >>> type(QQ_sympy().from_ZZ_sympy(Integer(2), ZZ_sympy()))
        <class 'sympy.core.numbers.Integer'>
        """
        return sympy_rat(a.p)

    def from_QQ_sympy(K1, a, K0):
        """
        Convert a SymPy `Rational` object to `dtype` (`Rational`).

        Example
        =======
        >>> from sympy import Rational
        >>> from sympy.polys.algebratools import QQ_sympy
        >>> QQ_sympy().from_QQ_sympy(Rational(3, 2), QQ_sympy())
        3/2
        >>> type(QQ_sympy().from_QQ_sympy(Rational(3, 2), QQ_sympy()))
        <class 'sympy.core.numbers.Rational'>
        """
        return a

    def from_ZZ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpz` object to `dtype` (`Rational`).

        Example
        =======
        >>> from sympy import Rational
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import QQ_sympy, ZZ_gmpy
        ...     a = QQ_sympy().from_ZZ_gmpy(gmpy.mpz(2), ZZ_gmpy())
        ... except ImportError:
        ...     a = Rational(2)
        >>> a
        2
        >>> type(a)
        <class 'sympy.core.numbers.Integer'>
        """
        return sympy_rat(python_int(a))

    def from_QQ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpq` object to `dtype` (`Rational`).

        Example
        =======
        >>> from sympy import Rational
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import QQ_sympy, QQ_gmpy
        ...     a = QQ_sympy().from_QQ_gmpy(gmpy.mpq(3, 2), QQ_gmpy())
        ... except ImportError:
        ...     a = Rational(3, 2)
        >>> a
        3/2
        >>> type(a)
        <class 'sympy.core.numbers.Rational'>
        """
        return sympy_rat(python_int(a.numer()),
                         python_int(a.denom()))

    def from_RR_sympy(K1, a, K0):
        """
        Convert a SymPy `Real` object to `dtype` (`Rational`).

        Example
        =======
        >>> from sympy import Real
        >>> from sympy.polys.algebratools import QQ_sympy, RR_sympy
        >>> QQ_sympy().from_RR_sympy(Real(1.5), RR_sympy())
        3/2
        >>> type(QQ_sympy().from_RR_sympy(Real(1.5), RR_sympy()))
        <class 'sympy.core.numbers.Rational'>
        """
        return sympy_rat(*K0.as_integer_ratio(a))

    def from_RR_mpmath(K1, a, K0):
        """
        Convert a mpmath `mpf` object to `dtype` (`Rational`).

        Example
        =======
        >>> from sympy import Rational
        >>> from sympy.mpmath import mpf
        >>> from sympy.polys.algebratools import QQ_sympy, RR_mpmath
        >>> QQ_sympy().from_RR_mpmath(mpf(1.5), RR_mpmath())
        3/2
        >>> type(QQ_sympy().from_RR_mpmath(mpf(1.5), RR_mpmath()))
        <class 'sympy.core.numbers.Rational'>
        """
        return sympy_rat(*K0.as_integer_ratio(a))

    def numer(self, a):
        """
        Returns the numerator of `a`.

        Example
        =======
        >>> from sympy import Rational
        >>> from sympy.polys.algebratools import QQ_sympy
        >>> QQ_sympy().numer(Rational(3, 2))
        3
        >>> type(QQ_sympy().numer(Rational(3, 2)))
        <class 'sympy.core.numbers.Integer'>
        """
        return sympy_int(a.p)

    def denom(self, a):
        """
        Returns the denominator of `a`.

        Example
        =======
        >>> from sympy import Rational
        >>> from sympy.polys.algebratools import QQ_sympy
        >>> QQ_sympy().denom(Rational(3, 2))
        2
        >>> type(QQ_sympy().denom(Rational(3, 2)))
        <class 'sympy.core.numbers.Integer'>
        """
        return sympy_int(a.q)

if HAS_GMPY:
    from gmpy import (
        mpz as gmpy_int,
        mpq as gmpy_rat,
        numer as gmpy_numer,
        denom as gmpy_denom,
        gcdext as gmpy_gcdex,
        gcd as gmpy_gcd,
        lcm as gmpy_lcm,
        sqrt as gmpy_sqrt,
        fac as gmpy_factorial,
    )

    class ZZ_gmpy(IntegerRing):
        """Integer ring based on GMPY mpz class. """

        dtype = gmpy_int
        zero  = dtype(0)
        one   = dtype(1)

        def __init__(self):
            pass

        def to_sympy(self, a):
            """
            Convert `a` to a SymPy object.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import ZZ_gmpy
            >>> ZZ_gmpy().to_sympy(gmpy.mpz(2))
            2
            >>> type(ZZ_gmpy().to_sympy(gmpy.mpz(2)))
            <class 'sympy.core.numbers.Integer'>
            """
            return sympy_int(int(a))

        def from_sympy(self, a):
            """
            Convert SymPy's Integer to `dtype` (`mpz`).

            Example
            =======
            >>> import gmpy
            >>> from sympy import S
            >>> from sympy.polys.algebratools import ZZ_gmpy
            >>> ZZ_gmpy().from_sympy(S.One)
            1
            >>> type(ZZ_gmpy().from_sympy(S.One))
            <type 'mpz'>
            """
            if a.is_Integer:
                return gmpy_int(a.p)
            elif a.is_Real and int(a) == a:
                return gmpy_int(int(a))
            else:
                raise CoercionFailed("expected Integer object, got %s" % a)

        def from_ZZ_python(K1, a, K0):
            """
            Convert a Python `int` object to `dtype` (`mpz`).

            Example
            =======
            >>> from sympy.polys.algebratools import ZZ_gmpy, ZZ_python
            >>> ZZ_gmpy().from_ZZ_python(2, ZZ_python())
            2
            >>> type(ZZ_gmpy().from_ZZ_python(2, ZZ_python()))
            <type 'mpz'>
            """
            return gmpy_int(a)

        def from_QQ_python(K1, a, K0):
            """
            Convert a Python `Fraction` object to `dtype` (`mpz`).

            Example
            =======
            >>> import gmpy
            >>> try:
            ...     import fractions # fractions only exists in >= Python 2.6
            ...     from sympy.polys.algebratools import ZZ_gmpy, QQ_python
            ...     a = ZZ_gmpy().from_QQ_python(fractions.Fraction(2, 1),
            ...     QQ_python())
            ... except ImportError:
            ...     a = gmpy.mpz(2)
            >>> a
            2
            >>> type(a)
            <type 'mpz'>
            """
            if a.denominator == 1:
                return gmpy_int(a.numerator)

        def from_ZZ_sympy(K1, a, K0):
            """
            Convert a SymPy `Integer` object to `dtype` (`mpz`).

            Example
            =======
            >>> from sympy import Integer
            >>> from sympy.polys.algebratools import ZZ_gmpy, ZZ_sympy
            >>> ZZ_gmpy().from_ZZ_sympy(Integer(2), ZZ_sympy())
            2
            >>> type(ZZ_gmpy().from_ZZ_sympy(Integer(2), ZZ_sympy()))
            <type 'mpz'>
            """
            return gmpy_int(a.p)

        def from_QQ_sympy(K1, a, K0):
            """
            Convert a SymPy `Rational` object to `dtype` (`mpz`).

            Example
            =======
            >>> from sympy import Rational
            >>> from sympy.polys.algebratools import ZZ_gmpy, QQ_sympy
            >>> ZZ_gmpy().from_QQ_sympy(Rational(2, 1), QQ_sympy())
            2
            >>> type(ZZ_gmpy().from_QQ_sympy(Rational(2, 1), QQ_sympy()))
            <type 'mpz'>
            """
            if a.q == 1:
                return gmpy_int(a.p)

        def from_ZZ_gmpy(K1, a, K0):
            """
            Convert a GMPY `mpz` object to `dtype` (`mpz`).

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import ZZ_gmpy
            >>> ZZ_gmpy().from_ZZ_gmpy(gmpy.mpz(2), ZZ_gmpy())
            2
            >>> type(ZZ_gmpy().from_ZZ_gmpy(gmpy.mpz(2), ZZ_gmpy()))
            <type 'mpz'>
            """
            return a

        def from_QQ_gmpy(K1, a, K0):
            """
            Convert a GMPY `mpq` object to `dtype` (`mpz`).

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import ZZ_gmpy, QQ_gmpy
            >>> ZZ_gmpy().from_QQ_gmpy(gmpy.mpq(2, 1), QQ_gmpy())
            2
            >>> type(ZZ_gmpy().from_QQ_gmpy(gmpy.mpq(2, 1), QQ_gmpy()))
            <type 'mpz'>
            """
            if a.denom() == 1:
                return a.numer()

        def from_RR_sympy(K1, a, K0):
            """
            Convert a SymPy `Real` object to `dtype` (`mpz`).

            Example
            =======
            >>> from sympy import Real
            >>> from sympy.polys.algebratools import ZZ_gmpy, RR_sympy
            >>> ZZ_gmpy().from_RR_sympy(Real(2.0), RR_sympy())
            2
            >>> type(ZZ_gmpy().from_RR_sympy(Real(2.0), RR_sympy()))
            <type 'mpz'>
            """
            p, q = K0.as_integer_ratio(a)

            if q == 1:
                return gmpy_int(p)

        def from_RR_mpmath(K1, a, K0):
            """
            Convert a mpmath `mpf` object to `dtype` (`mpz`).

            Example
            =======
            >>> from sympy.mpmath import mpf
            >>> from sympy.polys.algebratools import ZZ_gmpy, RR_mpmath
            >>> ZZ_gmpy().from_RR_mpmath(mpf(2.0), RR_mpmath())
            2
            >>> type(ZZ_gmpy().from_RR_mpmath(mpf(2.0), RR_mpmath()))
            <type 'mpz'>
            """
            p, q = K0.as_integer_ratio(a)

            if q == 1:
                return gmpy_int(p)

        def gcdex(self, a, b):
            """
            Extended GCD of `a` and `b`.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import ZZ_gmpy
            >>> ZZ_gmpy().gcdex(gmpy.mpz(12), gmpy.mpz(8))
            (1, -1, 4)
            >>> map(type, ZZ_gmpy().gcdex(gmpy.mpz(12), gmpy.mpz(8)))
            [<type 'mpz'>, <type 'mpz'>, <type 'mpz'>]
            """
            h, s, t = gmpy_gcdex(a, b)
            return s, t, h

        def gcd(self, a, b):
            """
            Returns the GCD of `a` and `b`.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import ZZ_gmpy
            >>> ZZ_gmpy().gcd(gmpy.mpz(12), gmpy.mpz(8))
            4
            >>> type(ZZ_gmpy().gcd(gmpy.mpz(12), gmpy.mpz(8)))
            <type 'mpz'>
            """
            return gmpy_gcd(a, b)

        def lcm(self, a, b):
            """
            Returns the LCM of `a` and `b`.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import ZZ_gmpy
            >>> ZZ_gmpy().lcm(gmpy.mpz(12), gmpy.mpz(8))
            24
            >>> type(ZZ_gmpy().lcm(gmpy.mpz(12), gmpy.mpz(8)))
            <type 'mpz'>
            """
            return gmpy_lcm(a, b)

        def sqrt(self, a):
            """
            Returns the square root of `a`.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import ZZ_gmpy
            >>> ZZ_gmpy().sqrt(gmpy.mpz(9))
            3
            >>> ZZ_gmpy().sqrt(gmpy.mpz(12))
            3
            >>> type(ZZ_gmpy().sqrt(gmpy.mpz(9)))
            <type 'mpz'>
            """
            return gmpy_sqrt(a)

        def factorial(self, a):
            """
            Returns the factorial of `a`.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import ZZ_gmpy
            >>> ZZ_gmpy().factorial(gmpy.mpz(4))
            24
            >>> type(ZZ_gmpy().factorial(gmpy.mpz(4)))
            <type 'mpz'>
            """
            return gmpy_factorial(int(a))

    class QQ_gmpy(RationalField):
        """Rational field based on GMPY mpq class. """

        dtype = gmpy_rat
        zero  = dtype(0)
        one   = dtype(1)

        def __init__(self):
            pass

        def to_sympy(self, a):
            """
            Convert `a` to a SymPy object.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import QQ_gmpy
            >>> QQ_gmpy().to_sympy(gmpy.mpq(3, 2))
            3/2
            >>> type(QQ_gmpy().to_sympy(gmpy.mpq(3, 2)))
            <class 'sympy.core.numbers.Rational'>
            """
            return sympy_rat(int(gmpy_numer(a)),
                             int(gmpy_denom(a)))

        def from_sympy(self, a):
            """
            Convert SymPy's Integer to `dtype` (`mpq`).

            Example
            =======
            >>> from sympy import Rational
            >>> import gmpy
            >>> from sympy.polys.algebratools import QQ_gmpy
            >>> QQ_gmpy().from_sympy(Rational(3, 2))
            3/2
            >>> type(QQ_gmpy().from_sympy(Rational(3, 2)))
            <type 'mpq'>
            """
            if a.is_Rational and a.q != 0:
                return gmpy_rat(a.p, a.q)
            elif a.is_Real:
                return gmpy_rat(*RR.as_integer_ratio(a))
            else:
                raise CoercionFailed("expected `Rational` object, got %s" % a)

        def from_ZZ_python(K1, a, K0):
            """
            Convert a Python `int` object to `dtype` (`mpq`).

            Example
            =======
            >>> from sympy.polys.algebratools import QQ_gmpy, ZZ_python
            >>> QQ_gmpy().from_ZZ_python(2, ZZ_python())
            2/1
            >>> type(QQ_gmpy().from_ZZ_python(2, ZZ_python()))
            <type 'mpq'>
            """
            return gmpy_rat(a)

        def from_QQ_python(K1, a, K0):
            """
            Convert a Python `Fraction` object to `dtype` (`mpq`).

            Example
            =======
            >>> import gmpy
            >>> try:
            ...     import fractions # fractions only exists in >= Python 2.6
            ...     from sympy.polys.algebratools import QQ_gmpy, QQ_python
            ...     a = QQ_gmpy().from_QQ_python(fractions.Fraction(3, 2),
            ...     QQ_python)
            ... except ImportError:
            ...     a = gmpy.mpq(3, 2)
            >>> a
            3/2
            >>> type(a)
            <type 'mpq'>
            """
            return gmpy_rat(a.numerator, a.denominator)

        def from_ZZ_sympy(K1, a, K0):
            """
            Convert a SymPy `Integer` object to `dtype` (`mpq`).

            Example
            =======
            >>> from sympy import Integer
            >>> from sympy.polys.algebratools import QQ_gmpy, ZZ_sympy
            >>> QQ_gmpy().from_ZZ_sympy(Integer(2), ZZ_sympy())
            2/1
            >>> type(QQ_gmpy().from_ZZ_sympy(Integer(2), ZZ_sympy()))
            <type 'mpq'>
            """
            return gmpy_rat(a.p)

        def from_QQ_sympy(K1, a, K0):
            """
            Convert a SymPy `Rational` object to `dtype` (`mpq`).

            Example
            =======
            >>> from sympy import Rational
            >>> from sympy.polys.algebratools import QQ_gmpy, QQ_sympy
            >>> QQ_gmpy().from_QQ_sympy(Rational(3, 2), QQ_sympy())
            3/2
            >>> type(QQ_gmpy().from_QQ_sympy(Rational(3, 2), QQ_sympy()))
            <type 'mpq'>
            """
            return gmpy_rat(a.p, a.q)

        def from_ZZ_gmpy(K1, a, K0):
            """
            Convert a GMPY `mpz` object to `dtype` (`mpq`).

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import QQ_gmpy, ZZ_gmpy
            >>> QQ_gmpy().from_ZZ_gmpy(gmpy.mpz(2), ZZ_gmpy())
            2/1
            >>> type(QQ_gmpy().from_ZZ_gmpy(gmpy.mpz(2), ZZ_gmpy()))
            <type 'mpq'>
            """
            return gmpy_rat(a)

        def from_QQ_gmpy(K1, a, K0):
            """
            Convert a GMPY `mpq` object to `dtype` (`mpq`).

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import QQ_gmpy
            >>> QQ_gmpy().from_QQ_gmpy(gmpy.mpq(3, 2), QQ_gmpy())
            3/2
            >>> type(QQ_gmpy().from_QQ_gmpy(gmpy.mpq(3, 2), QQ_gmpy()))
            <type 'mpq'>
            """
            return a

        def from_RR_sympy(K1, a, K0):
            """
            Convert a SymPy `Real` object to `dtype` (`mpq`).

            Example
            =======
            >>> from sympy import Real
            >>> from sympy.polys.algebratools import QQ_gmpy, RR_sympy
            >>> QQ_gmpy().from_RR_sympy(Real(1.5), RR_sympy())
            3/2
            >>> type(QQ_gmpy().from_RR_sympy(Real(1.5), RR_sympy()))
            <type 'mpq'>
            """
            return gmpy_rat(*K0.as_integer_ratio(a))

        def from_RR_mpmath(K1, a, K0):
            """
            Convert a mpmath `mpf` object to `dtype` (`mpq`).

            Example
            =======
            >>> from sympy.mpmath import mpf
            >>> from sympy.polys.algebratools import QQ_gmpy, RR_mpmath
            >>> QQ_gmpy().from_RR_mpmath(mpf(1.5), RR_mpmath())
            3/2
            >>> type(QQ_gmpy().from_RR_mpmath(mpf(1.5), RR_mpmath()))
            <type 'mpq'>
            """
            return gmpy_rat(*K0.as_integer_ratio(a))

        def exquo(self, a, b):
            """
            Exact quotient of `a` and `b`, implies `__div__`.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import QQ_gmpy
            >>> QQ_gmpy().exquo(gmpy.mpq(3, 2), gmpy.mpq(5, 6))
            9/5
            >>> type(QQ_gmpy().exquo(gmpy.mpq(3, 2), gmpy.mpq(5, 6)))
            <type 'mpq'>
            """
            return gmpy_rat(a.qdiv(b))

        def quo(self, a, b):
            """
            Quotient of `a` and `b`, implies `__div__`.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import QQ_gmpy
            >>> QQ_gmpy().quo(gmpy.mpq(3, 2), gmpy.mpq(5, 6))
            9/5
            >>> type(QQ_gmpy().quo(gmpy.mpq(3, 2), gmpy.mpq(5, 6)))
            <type 'mpq'>
            """
            return gmpy_rat(a.qdiv(b))

        def rem(self, a, b):
            """
            Remainder of `a` and `b`, implies nothing.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import QQ_gmpy
            >>> QQ_gmpy().rem(gmpy.mpq(3, 2), gmpy.mpq(5, 6))
            0/1
            >>> type(QQ_gmpy().rem(gmpy.mpq(3, 2), gmpy.mpq(5, 6)))
            <type 'mpq'>
            """
            return self.zero

        def div(self, a, b):
            """
            Division of `a` and `b`, implies `__div__`.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import QQ_gmpy
            >>> QQ_gmpy().div(gmpy.mpq(3, 2), gmpy.mpq(5, 6))
            (9/5, 0/1)
            >>> map(type, QQ_gmpy().div(gmpy.mpq(3, 2), gmpy.mpq(5, 6)))
            [<type 'mpq'>, <type 'mpq'>]
             """
            return gmpy_rat(a.qdiv(b)), self.zero

        def numer(self, a):
            """
            Returns the numerator of `a`.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import QQ_gmpy
            >>> QQ_gmpy().numer(gmpy.mpq(3, 2))
            3
            >>> type(QQ_gmpy().numer(gmpy.mpq(3, 2)))
            <type 'mpz'>
            """
            return gmpy_numer(a)

        def denom(self, a):
            """
            Returns the denominator of `a`.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import QQ_gmpy
            >>> QQ_gmpy().denom(gmpy.mpq(3, 2))
            2
            >>> type(QQ_gmpy().denom(gmpy.mpq(3, 2)))
            <type 'mpz'>
            """
            return gmpy_denom(a)

        def factorial(self, a):
            """
            Returns the factorial of `a`.

            Example
            =======
            >>> import gmpy
            >>> from sympy.polys.algebratools import QQ_gmpy
            >>> QQ_gmpy().factorial(gmpy.mpq(4, 1))
            24/1
            >>> type(QQ_gmpy().factorial(gmpy.mpq(4, 1)))
            <type 'mpq'>
            """
            return gmpy_rat(gmpy_factorial(int(a)))

class RealAlgebra(Algebra):
    """Abstract algebra for real numbers. """

    rep   = 'RR'

    is_Exact     = False
    is_Numerical = True

    has_CharacteristicZero = True

    def as_integer_ratio(self, a, **args):
        """
        Convert real number to a (numer, denom) pair.

        Example
        =======
        >>> from sympy.polys.algebratools import RR
        >>> RR.as_integer_ratio(1.5)
        (3, 2)
        """
        v, n = math.frexp(a) # XXX: hack, will work only for floats

        for i in xrange(300):
            if v != math.floor(v):
                v, n = 2*v, n-1
            else:
                break

        numer, denom = int(v), 1

        m = 1 << abs(n)

        if n > 0:
            numer *= m
        else:
            denom = m

        n, d = self.limit_denom(numer, denom, **args)

        if a and not n:
            return numer, denom
        else:
            return n, d

    def limit_denom(self, n, d, **args):
        """
        Find closest rational to `n/d` (up to `max_denom`).

        Example
        =======
        >>> from sympy.polys.algebratools import RR
        >>> RR.limit_denom(4521, 2045, max_denom=100)
        (42, 19)
        """
        max_denom = args.get('max_denom', 1000000)

        if d <= max_denom:
            return n, d

        self = QQ(n, d)

        p0, q0, p1, q1 = 0, 1, 1, 0

        while True:
            a  = n//d
            q2 = q0 + a*q1

            if q2 > max_denom:
                break

            p0, q0, p1, q1, n, d = \
                p1, q1, p0 + a*p1, q2, d, n - a*d

        k = (max_denom - q0)//q1

        P1, Q1 = p0 + k*p1, q0 + k*q1

        bound1 = QQ(P1, Q1)
        bound2 = QQ(p1, q1)

        if abs(bound2 - self) <= abs(bound1 - self):
            return p1, q1
        else:
            return P1, Q1

    def get_ring(self):
        """
        Returns a ring associated with `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import RR
        >>> RR.get_ring()
        Traceback (most recent call last):
        ...
        DomainError: there is no ring associated with RR
        """
        raise DomainError('there is no ring associated with %s' % self)

    def get_field(self):
        """
        Returns a field associated with `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import RR
        >>> RR.get_field()
        Traceback (most recent call last):
        ...
        DomainError: there is no field associated with RR
        """
        raise DomainError('there is no field associated with %s' % self)

    def get_exact(self):
        """
        Returns an exact domain associated with `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import RR
        >>> RR.get_exact()
        QQ
        """
        return QQ

    def exquo(self, a, b):
        """
        Exact quotient of `a` and `b`, implies `__div__`.

        Example
        =======
        >>> from sympy.polys.algebratools import RR
        >>> RR.exquo(RR(5.5), RR(5.0))
        1.1
        """
        return a / b

    def quo(self, a, b):
        """
        Quotient of `a` and `b`, implies `__div__`.

        Example
        =======
        >>> from sympy.polys.algebratools import RR
        >>> RR.quo(RR(5.5), RR(5.0))
        1.1
        """
        return a / b

    def rem(self, a, b):
        """
        Remainder of `a` and `b`, implies nothing.

        Example
        =======
        >>> from sympy.polys.algebratools import RR
        >>> RR.rem(RR(5.5), RR(5.0))
        0.0
        """
        return self.zero

    def div(self, a, b):
        """
        Division of `a` and `b`, implies `__div__`.

        Example
        =======
        >>> from sympy.polys.algebratools import RR
        >>> RR.div(RR(5.5), RR(5.0))
        (1.1, 0.0)
        """
        return a / b, self.zero

from sympy import (
    Real as sympy_mpf,
)

class RR_sympy(RealAlgebra):
    """Algebra for real numbers based on SymPy Real type. """

    dtype = sympy_mpf
    zero  = dtype(0)
    one   = dtype(1)

    def __init__(self):
        pass

    def to_sympy(self, a):
        """
        Convert `a` to a SymPy object.

        Example
        =======
        >>> from sympy import Real
        >>> from sympy.polys.algebratools import RR_sympy
        >>> RR_sympy().to_sympy(Real(1.5))
        1.50000000000000
        >>> type(RR_sympy().to_sympy(Real(1.5)))
        <class 'sympy.core.numbers.Real'>
        """
        return a

    def from_sympy(self, a):
        """
        Convert SymPy's Integer to `dtype` (`Real`).

        Example
        =======
        >>> from sympy import Real, Rational
        >>> from sympy.polys.algebratools import RR_sympy
        >>> RR_sympy().from_sympy(Rational(3, 2))
        1.50000000000000
        >>> type(RR_sympy().from_sympy(Rational(3, 2)))
        <class 'sympy.core.numbers.Real'>
        """
        b = a.evalf()

        if b.is_Real and b not in [S.Infinity, S.NegativeInfinity]:
            return b
        else:
            raise CoercionFailed("expected Real object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """
        Convert a Python `int` object to `dtype` (`Real`).

        Example
        =======
        >>> from sympy.polys.algebratools import RR_sympy, ZZ_python
        >>> RR_sympy().from_ZZ_python(2, ZZ_python())
        2
        >>> type(RR_sympy().from_ZZ_python(2, ZZ_python()))
        <class 'sympy.core.numbers.Integer'>
        """
        return sympy_mpf(a)

    def from_QQ_python(K1, a, K0):
        """
        Convert a Python `Fraction` object to `dtype` (`Real`).

        Example
        =======
        >>> from sympy import Rational
        >>> try:
        ...     import fractions # fractions only exists in >= Python 2.6
        ...     from sympy.polys.algebratools import RR_sympy, QQ_python
        ...     a = RR_sympy().from_QQ_python(fractions.Fraction(3, 2),
        ...     QQ_python())
        ... except ImportError:
        ...     a = Rational(3, 2)
        >>> a
        3/2
        >>> type(a)
        <class 'sympy.core.numbers.Rational'>
        """
        return sympy_mpf(a.numerator) / a.denominator

    def from_ZZ_sympy(K1, a, K0):
        """
        Convert a SymPy `Integer` object to `dtype` (`Real`).

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.polys.algebratools import RR_sympy, ZZ_sympy
        >>> RR_sympy().from_ZZ_sympy(Integer(2), ZZ_sympy())
        2
        >>> type(RR_sympy().from_ZZ_sympy(Integer(2), ZZ_sympy()))
        <class 'sympy.core.numbers.Integer'>
        """
        return sympy_mpf(a.p)

    def from_QQ_sympy(K1, a, K0):
        """
        Convert a SymPy `Rational` object to `dtype` (`Real`).

        Example
        =======
        >>> from sympy import Rational
        >>> from sympy.polys.algebratools import RR_sympy, QQ_sympy
        >>> RR_sympy().from_QQ_sympy(Rational(3, 2), QQ_sympy())
        1.50000000000000
        >>> type(RR_sympy().from_QQ_sympy(Rational(3, 2), QQ_sympy()))
        <class 'sympy.core.numbers.Real'>
        """
        return sympy_mpf(a)

    def from_ZZ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpz` object to `dtype` (`Real`).

        Example
        =======
        >>> from sympy import Real
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import RR_sympy, ZZ_gmpy
        ...     a = RR_sympy().from_ZZ_gmpy(gmpy.mpz(2), ZZ_gmpy())
        ... except ImportError:
        ...     a = Real(2)
        >>> a
        2
        >>> type(a)
        <class 'sympy.core.numbers.Integer'>
        """
        return sympy_mpf(int(a))

    def from_QQ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpq` object to `dtype` (`Real`).

        Example
        =======
        >>> from sympy import Rational
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import RR_sympy, QQ_gmpy
        ...     a = RR_sympy().from_QQ_gmpy(gmpy.mpq(3, 2), QQ_gmpy())
        ... except ImportError:
        ...     a = Rational(3, 2)
        >>> a
        3/2
        >>> type(a)
        <class 'sympy.core.numbers.Rational'>
        """
        return sympy_mpf(int(a.numer())) / int(a.denom())

    def from_RR_sympy(K1, a, K0):
        """
        Convert a SymPy `Real` object to `dtype` (`Real`).

        Example
        =======
        >>> from sympy import Real
        >>> from sympy.polys.algebratools import RR_sympy
        >>> RR_sympy().from_RR_sympy(Real(1.5), RR_sympy())
        1.50000000000000
        >>> type(RR_sympy().from_RR_sympy(Real(1.5), RR_sympy()))
        <class 'sympy.core.numbers.Real'>
        """
        return a

    def from_RR_mpmath(K1, a, K0):
        """
        Convert a mpmath `mpf` object to `dtype` (`Real`).

        Example
        =======
        >>> from sympy.mpmath import mpf
        >>> from sympy.polys.algebratools import RR_sympy, RR_mpmath
        >>> RR_sympy().from_RR_mpmath(mpf(1.5), RR_mpmath())
        1.50000000000000
        >>> type(RR_sympy().from_RR_mpmath(mpf(1.5), RR_mpmath()))
        <class 'sympy.core.numbers.Real'>
        """
        return sympy_mpf(a)

from sympy.mpmath import (
    mpf as mpmath_mpf,
)

class RR_mpmath(RealAlgebra):
    """Algebra for real numbers based on mpmath mpf type. """

    dtype = mpmath_mpf
    zero  = dtype(0)
    one   = dtype(1)

    def __init__(self):
        pass

    def to_sympy(self, a):
        """
        Convert `a` to a SymPy object.

        Example
        =======
        >>> from sympy.mpmath import mpf
        >>> from sympy.polys.algebratools import RR_mpmath
        >>> RR_mpmath().to_sympy(mpf(1.5))
        1.50000000000000
        >>> type(RR_mpmath().to_sympy(mpf(1.5)))
        <class 'sympy.core.numbers.Real'>
        """
        return sympy_mpf(a)

    def from_sympy(self, a):
        """
        Convert SymPy's Integer to `dtype` (`mpf`)

        Example
        =======
        >>> from sympy import Real
        >>> from sympy.polys.algebratools import RR_mpmath
        >>> RR_mpmath().from_sympy(Real(1.5))
        1.5
        >>> type(RR_mpmath().from_sympy(Real(1.5)))
        <class 'sympy.mpmath.ctx_mp_python.mpf'>
        """
        b = a.evalf()

        if b.is_Real and b not in [S.Infinity, S.NegativeInfinity]:
            return mpmath_mpf(b)
        else:
            raise CoercionFailed("expected Real object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """
        Convert a Python `int` object to `dtype` (`mpf`)

        Example
        =======
        >>> from sympy.polys.algebratools import RR_mpmath, ZZ_python
        >>> RR_mpmath().from_ZZ_python(2, ZZ_python())
        2.0
        >>> type(RR_mpmath().from_ZZ_python(2, ZZ_python()))
        <class 'sympy.mpmath.ctx_mp_python.mpf'>
        """
        return mpmath_mpf(a)

    def from_QQ_python(K1, a, K0):
        """
        Convert a Python `Fraction` object to `dtype` (`mpf`)

        Example
        =======
        >>> from sympy.mpmath import mpf
        >>> try:
        ...     import fractions # fractions only exists in >= Python 2.6
        ...     from sympy.polys.algebratools import RR_mpmath, QQ_python
        ...     a = RR_mpmath().from_QQ_python(fractions.Fraction(3, 2),
        ...     QQ_python())
        ... except ImportError:
        ...     a = mpf(1.5)
        >>> a
        1.5
        >>> type(a)
        <class 'sympy.mpmath.ctx_mp_python.mpf'>
        """
        return mpmath_mpf(a.numerator) / a.denominator

    def from_ZZ_sympy(K1, a, K0):
        """
        Convert a SymPy `Integer` object to `dtype` (`mpf`)

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.polys.algebratools import RR_mpmath, ZZ_sympy
        >>> RR_mpmath().from_ZZ_sympy(Integer(2), ZZ_sympy())
        2.0
        >>> type(RR_mpmath().from_ZZ_sympy(Integer(2), ZZ_sympy()))
        <class 'sympy.mpmath.ctx_mp_python.mpf'>
        """
        return mpmath_mpf(a.p)

    def from_QQ_sympy(K1, a, K0):
        """
        Convert a SymPy `Rational` object to `dtype` (`mpf`)

        Example
        =======
        >>> from sympy import Rational
        >>> from sympy.polys.algebratools import RR_mpmath, QQ_sympy
        >>> RR_mpmath().from_QQ_sympy(Rational(3, 2), QQ_sympy())
        1.5
        >>> type(RR_mpmath().from_QQ_sympy(Rational(3, 2), QQ_sympy()))
        <class 'sympy.mpmath.ctx_mp_python.mpf'>
        """
        return mpmath_mpf(a.p) / a.q

    def from_ZZ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpz` object to `dtype` (`mpf`)

        Example
        =======
        >>> from sympy.mpmath import mpf
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import RR_mpmath, ZZ_gmpy
        ...     a = RR_mpmath().from_ZZ_gmpy(gmpy.mpz(2), ZZ_gmpy())
        ... except ImportError:
        ...     a = mpf(2)
        >>> a
        2.0
        >>> type(a)
        <class 'sympy.mpmath.ctx_mp_python.mpf'>
        """
        return mpmath_mpf(int(a))

    def from_QQ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpq` object to `dtype` (`mpf`)

        Example
        =======
        >>> from sympy.mpmath import mpf
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import RR_mpmath, QQ_gmpy
        ...     a = RR_mpmath().from_QQ_gmpy(gmpy.mpq(3, 2), QQ_gmpy())
        ... except ImportError:
        ...     a = mpf(1.5)
        >>> a
        1.5
        >>> type(a)
        <class 'sympy.mpmath.ctx_mp_python.mpf'>
        """
        return mpmath_mpf(int(a.numer())) / int(a.denom())

    def from_RR_sympy(K1, a, K0):
        """
        Convert a SymPy `Real` object to `dtype` (`mpf`)

        Example
        =======
        >>> from sympy import Real
        >>> from sympy.polys.algebratools import RR_mpmath, RR_sympy
        >>> RR_mpmath().from_RR_sympy(Real(1.5), RR_sympy())
        1.5
        >>> type(RR_mpmath().from_RR_sympy(Real(1.5), RR_sympy()))
        <class 'sympy.mpmath.ctx_mp_python.mpf'>
        """
        return mpmath_mpf(a)

    def from_RR_mpmath(K1, a, K0):
        """
        Convert a mpmath `mpf` object to `dtype` (`mpf`)

        Example
        =======
        >>> from sympy.mpmath import mpf
        >>> from sympy.polys.algebratools import RR_mpmath
        >>> RR_mpmath().from_RR_mpmath(mpf(1.5), RR_mpmath())
        1.5
        >>> type(RR_mpmath().from_RR_mpmath(mpf(1.5), RR_mpmath()))
        <class 'sympy.mpmath.ctx_mp_python.mpf'>
        """
        return a

class FF_float(RealAlgebra): # XXX: tmp solution
    """Float domain. """

    rep   = 'FF'

    is_FF = True

    dtype = float
    zero  = dtype(0)
    one   = dtype(1)

    def __init__(self):
        pass

    def normal(self, a):
        """
        Normalize `a` with respect to `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import FF
        >>> FF.normal(1.0)
        1.0
        >>> FF.normal(1e-20)
        0.0
        """
        if abs(a) < 1e-15:
            return self.zero
        else:
            return self.dtype(a)

    def to_sympy(self, a):
        """
        Convert `a` to a SymPy object.

        Example
        =======
        >>> from sympy.polys.algebratools import FF
        >>> FF.to_sympy(1.5)
        1.50000000000000
        >>> type(FF.to_sympy(1.5))
        <class 'sympy.core.numbers.Real'>
        """
        return sympy_mpf(a)

    def from_sympy(self, a):
        """
        Convert SymPy's Integer to `dtype` (`float`).

        Example
        =======
        >>> from sympy import Real
        >>> from sympy.polys.algebratools import FF
        >>> FF.from_sympy(Real(1.5))
        1.5
        >>> type(FF.from_sympy(Real(1.5)))
        <type 'float'>
        """
        b = a.evalf()

        if b.is_Real and b not in [S.Infinity, S.NegativeInfinity]:
            return float(b)
        else:
            raise CoercionFailed("expected Real object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """
        Convert a Python `int` object to `dtype` (`float`).

        Example
        =======
        >>> from sympy.polys.algebratools import FF, ZZ_python
        >>> FF.from_ZZ_python(2, ZZ_python())
        2.0
        >>> type(FF.from_ZZ_python(2, ZZ_python()))
        <type 'float'>
        """
        return K1.dtype(a)

    def from_QQ_python(K1, a, K0):
        """
        Convert a Python `Fraction` object to `dtype` (`float`).

        Example
        =======
        >>> try:
        ...     import fractions # fractions only exists in >= Python 2.6
        ...     from sympy.polys.algebratools import FF, QQ_python
        ...     a = FF.from_QQ_python(fractions.Fraction(3, 2), QQ_python())
        ... except ImportError:
        ...     a = 1.5
        >>> a
        1.5
        >>> type(a)
        <type 'float'>
        """
        return K1.dtype(a.numerator) / a.denominator

    def from_ZZ_sympy(K1, a, K0):
        """
        Convert a SymPy `Integer` object to `dtype` (`float`).

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.polys.algebratools import FF, ZZ_sympy
        >>> FF.from_ZZ_sympy(Integer(2), ZZ_sympy())
        2.0
        >>> type(FF.from_ZZ_sympy(Integer(2), ZZ_sympy()))
        <type 'float'>
        """
        return K1.dtype(a.p)

    def from_QQ_sympy(K1, a, K0):
        """
        Convert a SymPy `Rational` object to `dtype` (`float`).

        Example
        =======
        >>> from sympy import Rational
        >>> from sympy.polys.algebratools import FF, QQ_sympy
        >>> FF.from_QQ_sympy(Rational(3, 2), QQ_sympy())
        1.5
        >>> type(FF.from_QQ_sympy(Rational(3, 2), QQ_sympy()))
        <type 'float'>
        """
        return K1.dtype(a.p) / a.q

    def from_ZZ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpz` object to `dtype` (`float`).

        Example
        =======
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import FF, ZZ_gmpy
        ...     a = FF.from_ZZ_gmpy(gmpy.mpz(2), ZZ_gmpy())
        ... except ImportError:
        ...     a = 2.0
        >>> a
        2.0
        >>> type(a)
        <type 'float'>
        """
        return K1.dtype(int(a))

    def from_QQ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpq` object to `dtype` (`float`).

        Example
        =======
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import FF, QQ_gmpy
        ...     a = FF.from_QQ_gmpy(gmpy.mpq(3, 2), QQ_gmpy())
        ... except ImportError:
        ...     a = 1.5
        >>> a
        1.5
        >>> type(a)
        <type 'float'>
        """
        return K1.dtype(int(a.numer())) / int(a.denom())

    def from_RR_sympy(K1, a, K0):
        """
        Convert a SymPy `Real` object to `dtype` (`float`).

        Example
        =======
        >>> from sympy import Real
        >>> from sympy.polys.algebratools import FF, RR_sympy
        >>> FF.from_RR_sympy(Real(1.5), RR_sympy())
        1.5
        >>> type(FF.from_RR_sympy(Real(1.5), RR_sympy()))
        <type 'float'>
        """
        return K1.dtype(a)

    def from_RR_mpmath(K1, a, K0):
        """
        Convert a mpmath `mpf` object to `dtype` (`float`).

        Example
        =======
        >>> from sympy.mpmath import mpf
        >>> from sympy.polys.algebratools import FF, RR_mpmath
        >>> FF.from_RR_mpmath(mpf(1.5), RR_mpmath())
        1.5
        >>> type(FF.from_RR_mpmath(mpf(1.5), RR_mpmath()))
        <type 'float'>
        """
        return K1.dtype(a)

    def complex_domain(self):
        """
        Returns a complex domain associated with `self`.

        Example
        =======
        >>> from sympy.polys.algebratools import FF
        >>> FF.complex_domain()
        CC
        """
        return CC

sympy_mpc = lambda a: sympy_mpf(a.real) + sympy_mpf(a.imag)*S.ImaginaryUnit

class CC_complex(RealAlgebra): # XXX: tmp solution
    """Complex domain. """

    rep   = 'CC'

    dtype = complex
    zero  = dtype(0)
    one   = dtype(1)

    def __init__(self):
        pass

    def to_sympy(self, a):
        """
        Convert `a` to a SymPy object.

        Example
        =======
        >>> from sympy.polys.algebratools import CC
        >>> CC.to_sympy(1+2j)
        1.0 + 2.0*I
        >>> type(CC.to_sympy(1+2j))
        <class 'sympy.core.add.Add'>
        """
        return sympy_mpc(a)

    def from_sympy(self, a):
        """
        Convert SymPy's Integer to `dtype` (`complex`)

        Example
        =======
        >>> from sympy import Real, I
        >>> from sympy.polys.algebratools import CC
        >>> CC.from_sympy(Real(1.0) + Real(2.0)*I)
        (1+2j)
        >>> type(CC.from_sympy(Real(1.0) + Real(2.0)*I))
        <type 'complex'>
        """
        b = a.evalf()
        r, i = a.as_real_imag()
        r, i = map(sympy_mpf, (r, i))
        if r.is_Real and r not in [S.Infinity, S.NegativeInfinity] or r is S.Zero:
            R = float(r)
        else:
            raise CoercionFailed("expected Real object, got %s" % a)

        if i.is_Real and i not in [S.Infinity, S.NegativeInfinity] or i is S.Zero:
            I = float(i)
        else:
            raise CoercionFailed("expected Real object, got %s" % a)

        return complex(R, I)

    def from_ZZ_python(K1, a, K0):
        """
        Convert a Python `int` object to `dtype` (`complex`)

        Example
        =======
        >>> from sympy.polys.algebratools import CC, ZZ_python
        >>> CC.from_ZZ_python(1, ZZ_python())
        (1+0j)
        >>> type(CC.from_ZZ_python(1, ZZ_python()))
        <type 'complex'>
        """
        return K1.dtype(a)

    def from_QQ_python(K1, a, K0):
        """
        Convert a Python `Fraction` object to `dtype` (`complex`)

        Example
        =======
        >>> try:
        ...     import fractions # fractions only exists in >= Python 2.6
        ...     from sympy.polys.algebratools import CC, QQ_python
        ...     a = CC.from_QQ_python(fractions.Fraction(3, 2), QQ_python())
        ... except ImportError:
        ...     a = 1.5+0j
        >>> a
        (1.5+0j)
        >>> type(a)
        <type 'complex'>
        """
        return K1.dtype(a.numerator) / a.denominator

    def from_ZZ_sympy(K1, a, K0):
        """
        Convert a SymPy `Integer` object to `dtype` (`complex`)

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.polys.algebratools import CC, ZZ_sympy
        >>> CC.from_ZZ_sympy(Integer(2), ZZ_sympy())
        (2+0j)
        >>> type(CC.from_ZZ_sympy(Integer(2), ZZ_sympy()))
        <type 'complex'>
        """
        return K1.dtype(a.p)

    def from_QQ_sympy(K1, a, K0):
        """
        Convert a SymPy `Rational` object to `dtype` (`complex`)

        Example
        =======
        >>> from sympy import Rational
        >>> from sympy.polys.algebratools import CC, QQ_sympy
        >>> CC.from_QQ_sympy(Rational(3, 2), QQ_sympy())
        (1.5+0j)
        >>> type(CC.from_QQ_sympy(Rational(3, 2), QQ_sympy()))
        <type 'complex'>
        """
        return K1.dtype(a.p) / a.q

    def from_ZZ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpz` object to `dtype` (`complex`)

        Example
        =======
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import CC, ZZ_gmpy
        ...     a = CC.from_ZZ_gmpy(gmpy.mpz(2), ZZ_gmpy())
        ... except ImportError:
        ...     a = 2+0j
        >>> a
        (2+0j)
        >>> type(a)
        <type 'complex'>
        """
        return K1.dtype(int(a))

    def from_QQ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpq` object to `dtype` (`complex`)

        Example
        =======
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import CC, QQ_gmpy
        ...     a = CC.from_QQ_gmpy(gmpy.mpq(3, 2), QQ_gmpy())
        ... except ImportError:
        ...     a = 1.5+0j
        >>> a
        (1.5+0j)
        >>> type(a)
        <type 'complex'>
        """
        return K1.dtype(int(a.numer())) / int(a.denom())

    def from_RR_sympy(K1, a, K0):
        """
        Convert a SymPy `Real` object to `dtype` (`complex`)

        Example
        =======
        >>> from sympy import Real
        >>> from sympy.polys.algebratools import CC, RR_sympy
        >>> CC.from_RR_sympy(Real(1.5), RR_sympy())
        (1.5+0j)
        >>> type(CC.from_RR_sympy(Real(1.5), RR_sympy()))
        <type 'complex'>
        """
        return K1.dtype(a)

    def from_RR_mpmath(K1, a, K0):
        """
        Convert a mpmath `mpf` object to `dtype` (`complex`)

        Example
        =======
        >>> from sympy.mpmath import mpf
        >>> from sympy.polys.algebratools import CC, RR_mpmath
        >>> CC.from_RR_mpmath(mpf(1.5), RR_mpmath())
        (1.5+0j)
        >>> type(CC.from_RR_mpmath(mpf(1.5), RR_mpmath()))
        <type 'complex'>
        """
        return K1.dtype(a)

    def from_FF_float(K1, a, K0):
        """
        Convert a `FF` object to `dtype` (`complex`)

        Example
        =======
        >>> from sympy.polys.algebratools import CC, FF
        >>> CC.from_FF_float(1.0, FF)
        (1+0j)
        >>> type(CC.from_FF_float(1.0, FF))
        <type 'complex'>
        """
        return K1.dtype(a)

    def real(self, a):
        """
        Returns the real part of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import CC
        >>> CC.real(CC(1+2j))
        1.0
        """
        return a.real

    def imag(self, a):
        """
        Returns the imaginary part of `a`.

        Example
        =======
        >>> from sympy.polys.algebratools import CC
        >>> CC.imag(CC(1+2j))
        2.0
        """
        return a.imag

FF = FF_float()
CC = CC_complex()

def _getenv(key, default=None):
    from os import getenv
    return getenv(key, default)

GROUND_TYPES = _getenv('SYMPY_GROUND_TYPES', 'gmpy').lower()

if GROUND_TYPES == 'python':  # XXX: needs 2.6 or better (at least for now)
    ZZ = ZZ_python()

    if HAS_FRACTION:
        QQ = QQ_python()
    elif HAS_GMPY:
        QQ = QQ_gmpy()
        GROUND_TYPES = 'python/gmpy'
    else:
        QQ = QQ_sympy()
        GROUND_TYPES = 'python/sympy'
elif GROUND_TYPES == 'sympy': # XXX: this is *very* slow, guess why ;)
    ZZ = ZZ_sympy()
    QQ = QQ_sympy()
elif GROUND_TYPES == 'gmpy':  # XXX: should be fine? sorry, but no, try -Qnew, damn
    if HAS_GMPY:
        ZZ = ZZ_gmpy()
        QQ = QQ_gmpy()
    else:
        ZZ = ZZ_python()

        if HAS_FRACTION:
            QQ = QQ_python()
            GROUND_TYPES = 'python'
        else:
            QQ = QQ_sympy()
            GROUND_TYPES = 'python/sympy'
else:
    raise ValueError("invalid ground types: %s" % GROUND_TYPES)

RR = RR_mpmath()

from sympy.polys.polyclasses import DMP, DMF, ANP

from sympy.polys.polyutils import (
    dict_from_basic,
    basic_from_dict,
    _dict_reorder,
)

class AlgebraicField(Field):
    """A class for representing algebraic number fields. """

    dtype        = ANP

    is_Numerical = True
    is_Algebraic = True

    has_assoc_Ring  = False
    has_assoc_Field = True

    has_CharacteristicZero = True

    def __init__(self, dom, *ext):
        if not dom.is_QQ:
            raise DomainError("ground domain must be a rational field")

        from sympy.polys.numberfields import to_number_field

        self.ext = to_number_field(ext)
        self.mod = self.ext.minpoly.rep

        self.dom  = dom

        self.gens = (self.ext,)
        self.unit = self([dom(1), dom(0)])

        self.zero = self.dtype.zero(self.mod.rep, dom)
        self.one  = self.dtype.one(self.mod.rep, dom)

    def __str__(self):
        return str(self.dom) + '<' + str(self.ext) + '>'

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.dom, self.ext))

    def __call__(self, a):
        """Construct an element of `self` domain from `a`. """
        return ANP(a, self.mod.rep, self.dom)

    def __eq__(self, other):
        """Returns `True` if two algebras are equivalent. """
        if self.dtype == other.dtype:
            return self.ext == other.ext
        else:
            return False

    def __ne__(self, other):
        """Returns `False` if two algebras are equivalent. """
        if self.dtype == other.dtype:
            return self.ext != other.ext
        else:
            return True

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        from sympy.polys.numberfields import AlgebraicNumber
        return AlgebraicNumber(self.ext, a).as_basic()

    def from_sympy(self, a):
        """Convert SymPy's expression to `dtype`. """
        try:
            return self([self.dom.from_sympy(a)])
        except CoercionFailed:
            pass

        from sympy.polys.numberfields import to_number_field

        try:
            return self(to_number_field(a, self.ext).native_coeffs())
        except (NotAlgebraic, IsomorphismFailed):
            raise CoercionFailed("%s is not a valid algebraic number in %s" % (a, self))

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return K1(K1.dom.convert(a, K0))

    def from_ZZ_sympy(K1, a, K0):
        """Convert a SymPy `Integer` object to `dtype`. """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_sympy(K1, a, K0):
        """Convert a SymPy `Rational` object to `dtype`. """
        return K1(K1.dom.convert(a, K0))

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return K1(K1.dom.convert(a, K0))

    def from_RR_sympy(K1, a, K0):
        """Convert a SymPy `Real` object to `dtype`. """
        return K1(K1.dom.convert(a, K0))

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        return K1(K1.dom.convert(a, K0))

    def get_ring(self):
        """Returns a ring associated with `self`. """
        raise DomainError('there is no ring associated with %s' % self)

    def is_positive(self, a):
        """Returns True if `a` is positive. """
        return self.dom.is_positive(a.LC())

    def is_negative(self, a):
        """Returns True if `a` is negative. """
        return self.dom.is_negative(a.LC())

    def is_nonpositive(self, a):
        """Returns True if `a` is non-positive. """
        return self.dom.is_nonpositive(a.LC())

    def is_nonnegative(self, a):
        """Returns True if `a` is non-negative. """
        return self.dom.is_nonnegative(a.LC())

    def numer(self, a):
        """Returns the numerator of `a`. """
        return a

    def denom(self, a):
        """Returns the denominator of `a`. """
        return self.one

class PolynomialRing(Ring):
    """A class for representing multivariate polynomial rings. """

    dtype        = DMP
    is_Poly      = True
    is_Composite = True

    has_assoc_Ring         = True
    has_assoc_Field        = True

    has_CharacteristicZero = True

    def __init__(self, dom, *gens):
        if not gens:
            raise GeneratorsNeeded("generators not specified")

        lev = len(gens) - 1

        self.zero = self.dtype.zero(lev, dom)
        self.one  = self.dtype.one(lev, dom)

        self.dom  = dom
        self.gens = gens

    def __str__(self):
        return str(self.dom) + '[' + ','.join(map(str, self.gens)) + ']'

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.dom, self.gens))

    def __call__(self, a):
        """Construct an element of `self` domain from `a`. """
        return DMP(a, self.dom, len(self.gens)-1)

    def __eq__(self, other):
        """Returns `True` if two algebras are equivalent. """
        if self.dtype == other.dtype:
            return self.gens == other.gens and self.dom == other.dom
        else:
            return False

    def __ne__(self, other):
        """Returns `False` if two algebras are equivalent. """
        if self.dtype == other.dtype:
            return self.gens != other.gens
        else:
            return True

    def to_sympy(self, a):
        """
        Convert `a` to a SymPy object.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ[x].to_sympy(DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ))
        1 + x**2
        >>> type(ZZ[x].to_sympy(DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ)))
        <class 'sympy.core.add.Add'>
        """
        return basic_from_dict(a.to_sympy_dict(), *self.gens)

    def from_sympy(self, a):
        r"""
        Convert SymPy's expression to `dtype` (`DMP`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ[x].from_sympy(x**2 + 1) == \
        ... DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ)
        True
        >>> type(ZZ[x].from_sympy(x**2 + 1))
        <class 'sympy.polys.polyclasses.DMP'>
        """
        rep = dict_from_basic(a, self.gens)

        for k, v in rep.iteritems():
            rep[k] = self.dom.from_sympy(v)

        return self(rep)

    def from_ZZ_python(K1, a, K0):
        r"""
        Convert a Python `int` object to `dtype` (`DMP`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ, ZZ_python
        >>> ZZ[x].from_ZZ_python(1, ZZ_python()) == \
        ... DMP([ZZ(1)], ZZ)
        True
        >>> type(ZZ[x].from_ZZ_python(1, ZZ_python()))
        <class 'sympy.polys.polyclasses.DMP'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_python(K1, a, K0):
        r"""
        Convert a Python `Fraction` object to `dtype` (`DMP`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> try:
        ...     import fractions # fractions only exists in >= Python 2.6
        ...     from sympy.polys.algebratools import QQ_python
        ...     a = QQ[x].from_QQ_python(fractions.Fraction(3, 2), QQ_python())
        ... except ImportError:
        ...     a = DMP([QQ(3, 2)], QQ)
        >>> a == \
        ... DMP([QQ(3, 2)], QQ)
        True
        >>> type(a)
        <class 'sympy.polys.polyclasses.DMP'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_ZZ_sympy(K1, a, K0):
        r"""
        Convert a SymPy `Integer` object to `dtype` (`DMP`).

        Example
        =======
        >>> from sympy import Integer
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ, ZZ_sympy
        >>> ZZ[x].from_ZZ_sympy(Integer(2), ZZ_sympy()) == \
        ... DMP([ZZ(2)], ZZ)
        True
        >>> type(ZZ[x].from_ZZ_sympy(Integer(2), ZZ_sympy()))
        <class 'sympy.polys.polyclasses.DMP'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_sympy(K1, a, K0):
        r"""
        Convert a SymPy `Rational` object to `dtype` (`DMP`).

        Example
        =======
        >>> from sympy import Rational
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ, QQ_sympy
        >>> QQ[x].from_QQ_sympy(Rational(3, 2), QQ_sympy()) == \
        ... DMP([QQ(3, 2)], QQ)
        True
        >>> type(QQ[x].from_QQ_sympy(Rational(3, 2), QQ_sympy()))
        <class 'sympy.polys.polyclasses.DMP'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_ZZ_gmpy(K1, a, K0):
        r"""
        Convert a GMPY `mpz` object to `dtype` (`DMP`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import ZZ_gmpy
        ...     a = ZZ[x].from_ZZ_gmpy(gmpy.mpz(2), ZZ_gmpy())
        ... except ImportError:
        ...     a = DMP([ZZ(2)], ZZ)
        >>> a == \
        ... DMP([ZZ(2)], ZZ)
        True
        >>> type(a)
        <class 'sympy.polys.polyclasses.DMP'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_gmpy(K1, a, K0):
        r"""
        Convert a GMPY `mpq` object to `dtype` (`DMP`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import QQ_gmpy
        ...     a = QQ[x].from_QQ_gmpy(gmpy.mpq(3, 2), QQ_gmpy())
        ... except ImportError:
        ...     a = DMP([QQ(3, 2)], QQ)
        >>> a == \
        ... DMP([QQ(3, 2)], QQ)
        True
        >>> type(a)
        <class 'sympy.polys.polyclasses.DMP'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_RR_sympy(K1, a, K0):
        r"""
        Convert a SymPy `Real` object to `dtype` (`DMP`).

        Example
        =======
        >>> from sympy import Real
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ, RR_sympy
        >>> QQ[x].from_RR_sympy(Real(1.5), RR_sympy()) == \
        ... DMP([QQ(3, 2)], QQ)
        True
        >>> type(QQ[x].from_RR_sympy(Real(1.5), RR_sympy()))
        <class 'sympy.polys.polyclasses.DMP'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_RR_mpmath(K1, a, K0):
        r"""
        Convert a mpmath `mpf` object to `dtype` (`DMP`).

        Example
        =======
        >>> from sympy.mpmath import mpf
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ, RR_mpmath
        >>> QQ[x].from_RR_mpmath(mpf(1.5), RR_mpmath()) == \
        ... DMP([QQ(3, 2)], QQ)
        True
        >>> type(QQ[x].from_RR_mpmath(mpf(1.5), RR_mpmath()))
        <class 'sympy.polys.polyclasses.DMP'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_PolynomialRing(K1, a, K0):
        """
        Convert a `DMP` object to `dtype` (`DMP`).

        Example
        =======
        >>> from sympy.abc import x, y
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> a = ZZ[x].from_PolynomialRing(DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ[y]),
        ... ZZ[y])
        >>> a == DMP([1], ZZ)
        True
        >>> b = ZZ[x].from_PolynomialRing(DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ[x]),
        ... ZZ[x])
        >>> b == DMP([1, 0, 1], ZZ[x])
        True
        >>> (type(a), type(b))
        (<class 'sympy.polys.polyclasses.DMP'>,
         <class 'sympy.polys.polyclasses.DMP'>)
        """
        if K1.gens == K0.gens:
            if K1.dom == K0.dom:
                return a
            else:
                return a.convert(K1.dom)
        else:
            monoms, coeffs = _dict_reorder(a.to_dict(), K0.gens, K1.gens)

            if K1.dom != K0.dom:
                coeffs = [ K1.dom.convert(c, K0.dom) for c in coeffs ]

            return K1(dict(zip(monoms, coeffs)))

    def from_FractionField(K1, a, K0):
        """
        Convert a `DMF` object to `dtype` (`DMP`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP, DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> a = ZZ[x].from_FractionField(DMF(([ZZ(1), ZZ(1)], [ZZ(1)]), ZZ),
        ... ZZ[x])
        >>> a == DMP([ZZ(1), ZZ(1)], ZZ)
        True
        >>> type(a)
        <class 'sympy.polys.polyclasses.DMP'>
        """
        if a.denom().is_one:
            return K1.from_PolynomialRing(a.numer(), K0)

    def get_field(self):
        """
        Returns a field associated with `self`.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ[x].get_field()
        ZZ(x)
        >>> type(ZZ[x].get_field())
        <class 'sympy.polys.algebratools.FractionField'>
        """
        return FractionField(self.dom, *self.gens)

    def poly_ring(self, *gens):
        """
        Returns a polynomial ring, i.e. `K[X]`.

        Example
        =======
        >>> from sympy.abc import x, y
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ[x].poly_ring(y)
        Traceback (most recent call last):
        ...
        NotImplementedError: nested domains not allowed
        """
        raise NotImplementedError('nested domains not allowed')

    def frac_field(self, *gens):
        """
        Returns a fraction field, i.e. `K(X)`.

        Example
        =======
        >>> from sympy.abc import x, y
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ[x].frac_field(y)
        Traceback (most recent call last):
        ...
        NotImplementedError: nested domains not allowed
        """
        raise NotImplementedError('nested domains not allowed')

    def is_positive(self, a):
        """
        Returns True if `LC(a)` is positive.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ[x].is_positive(DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ))
        True
        >>> ZZ[x].is_positive(DMP([ZZ(-1), ZZ(0), ZZ(1)], ZZ))
        False
        >>> ZZ[x].is_positive(DMP([], ZZ))
        False
        """
        return self.dom.is_positive(a.LC())

    def is_negative(self, a):
        """
        Returns True if `LC(a)` is negative.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ[x].is_negative(DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ))
        False
        >>> ZZ[x].is_negative(DMP([ZZ(-1), ZZ(0), ZZ(1)], ZZ))
        True
        >>> ZZ[x].is_negative(DMP([], ZZ))
        False
        """
        return self.dom.is_negative(a.LC())

    def is_nonpositive(self, a):
        """
        Returns True if `LC(a)` is non-positive.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ[x].is_nonpositive(DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ))
        False
        >>> ZZ[x].is_nonpositive(DMP([ZZ(-1), ZZ(0), ZZ(1)], ZZ))
        True
        >>> ZZ[x].is_nonpositive(DMP([], ZZ))
        True
        """
        return self.dom.is_nonpositive(a.LC())

    def is_nonnegative(self, a):
        """
        Returns True if `LC(a)` is non-negative.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ[x].is_nonnegative(DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ))
        True
        >>> ZZ[x].is_nonnegative(DMP([ZZ(-1), ZZ(0), ZZ(1)], ZZ))
        False
        >>> ZZ[x].is_nonnegative(DMP([], ZZ))
        True
        """
        return self.dom.is_nonnegative(a.LC())

    def gcdex(self, a, b):
        r"""
        Extended GCD of `a` and `b`.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> QQ[x].gcdex(DMP([QQ(1), QQ(0), QQ(-1)], QQ),
        ... DMP([QQ(1), QQ(1), QQ(0)], QQ)) == \
        ... (DMP([QQ(-1)], QQ), DMP([QQ(1)], QQ), DMP([QQ(1), QQ(1)], QQ))
        True
        """
        return a.gcdex(b)

    def gcd(self, a, b):
        r"""
        Returns the GCD of `a` and `b`.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> QQ[x].gcd(DMP([QQ(1), QQ(0), QQ(-1)], QQ),
        ... DMP([QQ(1), QQ(1), QQ(0)], QQ)) == \
        ... DMP([QQ(1), QQ(1)], QQ)
        True
        """
        return a.gcd(b)

    def lcm(self, a, b):
        r"""
        Returns the LCM of `a` and `b`.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> QQ[x].lcm(DMP([QQ(1), QQ(0), QQ(-1)], QQ),
        ... DMP([QQ(1), QQ(1), QQ(0)], QQ)) == \
        ... DMP([QQ(1), QQ(0), QQ(-1), QQ(0)], QQ)
        True
        """
        return a.lcm(b)

    def factorial(self, a):
        """Returns the factorial of `a`. """
        # XXX: What is this supposed to do?
        return self.dtype(self.dom.factorial(a))

class FractionField(Field):
    """A class for representing rational function fields. """

    dtype        = DMF
    is_Frac      = True
    is_Composite = True

    has_assoc_Ring         = True
    has_assoc_Field        = True

    has_CharacteristicZero = True

    def __init__(self, dom, *gens):
        if not gens:
            raise GeneratorsNeeded("generators not specified")

        lev = len(gens) - 1

        self.zero = self.dtype.zero(lev, dom)
        self.one  = self.dtype.one(lev, dom)

        self.dom  = dom
        self.gens = gens

    def __str__(self):
        return str(self.dom) + '(' + ','.join(map(str, self.gens)) + ')'

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.dom, self.gens))

    def __call__(self, a):
        """Construct an element of `self` domain from `a`. """
        return DMF(a, self.dom, len(self.gens)-1)

    def __eq__(self, other):
        """Returns `True` if two algebras are equivalent. """
        if self.dtype == other.dtype:
            return self.gens == other.gens and self.dom == other.dom
        else:
            return False

    def __ne__(self, other):
        """Returns `False` if two algebras are equivalent. """
        if self.dtype == other.dtype:
            return self.gens != other.gens
        else:
            return True

    def to_sympy(self, a):
        """
        Convert `a` to a SymPy object.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.frac_field(x).to_sympy(DMF(([ZZ(1), ZZ(2)], [ZZ(1), ZZ(1)]), ZZ))
        (2 + x)/(1 + x)
        >>> type(ZZ.frac_field(x).to_sympy(DMF(([ZZ(1), ZZ(2)],
        ... [ZZ(1), ZZ(1)]), ZZ)))
        <class 'sympy.core.mul.Mul'>
        """
        return (basic_from_dict(a.numer().to_sympy_dict(), *self.gens) /
                basic_from_dict(a.denom().to_sympy_dict(), *self.gens))

    def from_sympy(self, a):
        r"""
        Convert SymPy's expression to `dtype` (`DMF`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.frac_field(x).from_sympy((2 + x)/(1 + x)) == \
        ... DMF(([ZZ(1), ZZ(2)], [ZZ(1), ZZ(1)]), ZZ)
        True
        >>> type(ZZ.frac_field(x).from_sympy((2 + x)/(2 + x)))
        <class 'sympy.polys.polyclasses.DMF'>
        """
        p, q = a.as_numer_denom()

        num = dict_from_basic(p, self.gens)
        den = dict_from_basic(q, self.gens)

        for k, v in num.iteritems():
            num[k] = self.dom.from_sympy(v)

        for k, v in den.iteritems():
            den[k] = self.dom.from_sympy(v)

        return self((num, den)).cancel()

    def from_ZZ_python(K1, a, K0):
        r"""
        Convert a Python `int` object to `dtype` (`DMF`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ, ZZ_python
        >>> ZZ.frac_field(x).from_ZZ_python(2, ZZ_python()) == \
        ... DMF(([ZZ(2)], [ZZ(1)]), ZZ)
        True
        >>> type(ZZ.frac_field(x).from_ZZ_python(2, ZZ_python()))
        <class 'sympy.polys.polyclasses.DMF'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_python(K1, a, K0):
        """
        Convert a Python `Fraction` object to `dtype` (`DMF`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import QQ
        >>> try:
        ...     import fractions # fractions only exists in >= Python 2.6
        ...     from sympy.polys.algebratools import QQ_python
        ...     a = QQ.frac_field(x).from_QQ_python(fractions.Fraction(3, 2),
        ...     QQ_python())
        ... except ImportError:
        ...     a = DMF(([QQ(3, 2)], [QQ(1)]), QQ)
        >>> a == DMF(([QQ(3, 2)], [QQ(1)]), QQ)
        True
        >>> type(a)
        <class 'sympy.polys.polyclasses.DMF'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_ZZ_sympy(K1, a, K0):
        r"""
        Convert a SymPy `Integer` object to `dtype` (`DMF`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy import Integer
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ, ZZ_sympy
        >>> ZZ.frac_field(x).from_ZZ_sympy(Integer(2), ZZ_sympy()) == \
        ... DMF(([ZZ(2)], [ZZ(1)]), ZZ)
        True
        >>> type(ZZ.frac_field(x).from_ZZ_sympy(Integer(2), ZZ_sympy()))
        <class 'sympy.polys.polyclasses.DMF'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_sympy(K1, a, K0):
        r"""
        Convert a SymPy `Rational` object to `dtype` (`DMF`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy import Rational
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import QQ, QQ_sympy
        >>> QQ.frac_field(x).from_QQ_sympy(Rational(3, 2), QQ_sympy()) == \
        ... DMF(([QQ(3, 2)], [QQ(1)]), QQ)
        True
        >>> type(QQ.frac_field(x).from_QQ_sympy(Rational(3, 2), QQ_sympy()))
        <class 'sympy.polys.polyclasses.DMF'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_ZZ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpz` object to `dtype` (`DMF`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import ZZ_gmpy
        ...     a = ZZ.frac_field(x).from_ZZ_gmpy(gmpy.mpz(2), ZZ_gmpy())
        ... except ImportError:
        ...     a = DMF(([ZZ(2)], [ZZ(1)]), ZZ)
        >>> a == DMF(([ZZ(2)], [ZZ(1)]), ZZ)
        True
        >>> type(a)
        <class 'sympy.polys.polyclasses.DMF'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_gmpy(K1, a, K0):
        """
        Convert a GMPY `mpq` object to `dtype` (`DMF`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import QQ
        >>> try:
        ...     import gmpy
        ...     from sympy.polys.algebratools import QQ_gmpy
        ...     a = QQ.frac_field(x).from_QQ_gmpy(gmpy.mpq(3, 2), QQ_gmpy())
        ... except ImportError:
        ...     a = DMF(([QQ(3, 2)], [QQ(1)]), QQ)
        >>> a == DMF(([QQ(3, 2)], [QQ(1)]), QQ)
        True
        >>> type(a)
        <class 'sympy.polys.polyclasses.DMF'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_RR_sympy(K1, a, K0):
        r"""
        Convert a SymPy `Real` object to `dtype`(`DMF`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy import Real
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import QQ, RR_sympy
        >>> QQ.frac_field(x).from_RR_sympy(Real(1.5), RR_sympy()) == \
        ... DMF(([QQ(3, 2)], [QQ(1)]), QQ)
        True
        >>> type(QQ.frac_field(x).from_RR_sympy(Real(1.5), RR_sympy()))
        <class 'sympy.polys.polyclasses.DMF'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_RR_mpmath(K1, a, K0):
        r"""
        Convert a mpmath `mpf` object to `dtype`(`DMF`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.mpmath import mpf
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import QQ, RR_mpmath
        >>> QQ.frac_field(x).from_RR_mpmath(mpf(1.5), RR_mpmath()) == \
        ... DMF(([QQ(3, 2)], [QQ(1)]), QQ)
        True
        >>> type(QQ.frac_field(x).from_RR_mpmath(mpf(1.5), RR_mpmath()))
        <class 'sympy.polys.polyclasses.DMF'>
        """
        return K1(K1.dom.convert(a, K0))

    def from_PolynomialRing(K1, a, K0):
        r"""
        Convert a `DMP` object to `dtype`(`DMF`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF, DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.frac_field(x).from_PolynomialRing(DMP([ZZ(1), ZZ(-1)], ZZ),
        ... ZZ[x]) == \
        ... DMF(([ZZ(1), ZZ(-1)], [ZZ(1)]), ZZ)
        True
        >>> type(ZZ.frac_field(x).from_PolynomialRing(DMP([ZZ(1), ZZ(-1)],
        ... ZZ), ZZ[x]))
        <class 'sympy.polys.polyclasses.DMF'>
        """
        if K1.gens == K0.gens:
            if K1.dom == K0.dom:
                return K1(a.rep)
            else:
                return K1(a.convert(K1.dom).rep)
        else:
            monoms, coeffs = _dict_reorder(a.to_dict(), K0.gens, K1.gens)

            if K1.dom != K0.dom:
                coeffs = [ K1.dom.convert(c, K0.dom) for c in coeffs ]

            return K1(dict(zip(monoms, coeffs)))

    def from_FractionField(K1, a, K0):
        r"""
        Convert a `DMF` object to `dtype`(`DMF`).

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.frac_field(x).from_FractionField(DMF(([ZZ(1), ZZ(2)],
        ... [ZZ(1), ZZ(1)]), ZZ), ZZ[x]) == \
        ... DMF(([ZZ(1), ZZ(2)], [ZZ(1), ZZ(1)]), ZZ)
        True
        >>> type(ZZ.frac_field(x).from_FractionField(DMF(([ZZ(1), ZZ(2)],
        ... [ZZ(1), ZZ(2)]), ZZ), ZZ[x]))
        <class 'sympy.polys.polyclasses.DMF'>
        """
        if K1.gens == K0.gens:
            if K1.dom == K0.dom:
                return a
            else:
                return K1(a.numer().convert(K1.dom).rep,
                          a.denom().convert(K1.dom).rep)

    def get_ring(self):
        """
        Returns a ring associated with `self`.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.frac_field(x).get_ring()
        ZZ[x]
        """
        return PolynomialRing(self.dom, *self.gens)

    def poly_ring(self, *gens):
        """
        Returns a polynomial ring, i.e. `K[X]`.

        Example
        =======
        >>> from sympy.abc import x, y
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.frac_field(x).poly_ring(y)
        Traceback (most recent call last):
        ...
        NotImplementedError: nested domains not allowed
        """
        raise NotImplementedError('nested domains not allowed')

    def frac_field(self, *gens):
        """
        Returns a fraction field, i.e. `K(X)`.

        Example
        =======
        >>> from sympy.abc import x, y
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.frac_field(x).frac_field(y)
        Traceback (most recent call last):
        ...
        NotImplementedError: nested domains not allowed
        """
        raise NotImplementedError('nested domains not allowed')

    def is_positive(self, a):
        """
        Returns True if `a` is positive.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.frac_field(x).is_positive(DMF(([ZZ(1), ZZ(2)],
        ... [ZZ(1), ZZ(2)]), ZZ))
        True
        >>> ZZ.frac_field(x).is_positive(DMF(([ZZ(-1), ZZ(2)],
        ... [ZZ(1), ZZ(2)]), ZZ))
        False
        >>> ZZ.frac_field(x).is_positive(DMF(([], [ZZ(1)]), ZZ))
        False
        """
        return self.dom.is_positive(a.numer().LC())

    def is_negative(self, a):
        """
        Returns True if `a` is negative.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.frac_field(x).is_negative(DMF(([ZZ(1), ZZ(2)],
        ... [ZZ(1), ZZ(2)]), ZZ))
        False
        >>> ZZ.frac_field(x).is_negative(DMF(([ZZ(-1), ZZ(2)],
        ... [ZZ(1), ZZ(2)]), ZZ))
        True
        >>> ZZ.frac_field(x).is_negative(DMF(([], [ZZ(1)]), ZZ))
        False
        """
        return self.dom.is_negative(a.numer().LC())

    def is_nonpositive(self, a):
        """
        Returns True if `a` is non-positive.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.frac_field(x).is_nonpositive(DMF(([ZZ(1), ZZ(2)],
        ... [ZZ(1), ZZ(2)]), ZZ))
        False
        >>> ZZ.frac_field(x).is_nonpositive(DMF(([ZZ(-1), ZZ(2)],
        ... [ZZ(1), ZZ(2)]), ZZ))
        True
        >>> ZZ.frac_field(x).is_nonpositive(DMF(([], [ZZ(1)]), ZZ))
        True
        """
        return self.dom.is_nonpositive(a.numer().LC())

    def is_nonnegative(self, a):
        """
        Returns True if `a` is non-negative.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.frac_field(x).is_nonnegative(DMF(([ZZ(1), ZZ(2)],
        ... [ZZ(1), ZZ(2)]), ZZ))
        True
        >>> ZZ.frac_field(x).is_nonnegative(DMF(([ZZ(-1), ZZ(2)],
        ... [ZZ(1), ZZ(2)]), ZZ))
        False
        >>> ZZ.frac_field(x).is_nonnegative(DMF(([], [ZZ(1)]), ZZ))
        True
        """
        return self.dom.is_nonnegative(a.numer().LC())

    def numer(self, a):
        r"""
        Returns the numerator of `a`.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF, DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.frac_field(x).numer(DMF(([ZZ(1), ZZ(2)],
        ... [ZZ(1), ZZ(1)]), ZZ)) == \
        ... DMP([ZZ(1), ZZ(2)], ZZ)
        True
        """
        return a.numer()

    def denom(self, a):
        r"""
        Returns the denominator of `a`.

        Example
        =======
        >>> from sympy.abc import x
        >>> from sympy.polys.polyclasses import DMF, DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> ZZ.frac_field(x).denom(DMF(([ZZ(1), ZZ(2)],
        ... [ZZ(1), ZZ(1)]), ZZ)) == \
        ... DMP([ZZ(1), ZZ(1)], ZZ)
        True
        """
        return a.denom()

    def factorial(self, a):
        """Returns the factorial of `a`. """
        # XXX: What is this suppposed to do?
        return self.dtype(self.dom.factorial(a))

class ExpressionDomain(Field):
    """A class for arbitrary expressions. """

    is_EX = True

    class Expression(object):
        """An arbitrary expression. """

        __slots__ = ['ex']

        def __init__(self, ex):
            if not isinstance(ex, self.__class__):
                self.ex = sympify(ex)
            else:
                self.ex = ex.ex

        def __repr__(self):
            return 'EX(%s)' % repr(self.ex)

        def __str__(self):
            return 'EX(%s)' %str(self.ex)

        def __hash__(self):
            return hash((self.__class__.__name__, self.ex))

        def as_basic(self):
            return self.ex

        def numer(self):
            return EX(self.ex.as_numer_denom()[0])

        def denom(self):
            return EX(self.ex.as_numer_denom()[1])

        def simplify(self, ex):
            from sympy import cancel
            return self.__class__(cancel(ex))

        def __abs__(self):
            return self.__class__(abs(self.ex))

        def __neg__(self):
            return self.__class__(-self.ex)

        def __add__(self, g):
            return self.simplify(self.ex+self.__class__(g).ex)

        def __radd__(self, g):
            return self.simplify(self.__class__(g).ex+self.ex)

        def __sub__(self, g):
            return self.simplify(self.ex-self.__class__(g).ex)

        def __rsub__(self, g):
            return self.simplify(self.__class__(g).ex-self.ex)

        def __mul__(self, g):
            return self.simplify(self.ex*self.__class__(g).ex)

        def __rmul__(self, g):
            return self.simplify(self.__class__(g).ex*self.ex)

        def __pow__(self, n):
            return self.simplify(self.ex**n)

        def __div__(self, g):
            return self.simplify(self.ex/self.__class__(g).ex)

        def __rdiv__(self, g):
            return self.simplify(self.__class__(g).ex/self.ex)

        def __truediv__(self, g):
            return self.simplify(self.ex/self.__class__(g).ex)

        def __rtruediv__(self, g):
            return self.simplify(self.__class__(g).ex/self.ex)

        def __eq__(self, g):
            return self.ex == self.__class__(g).ex

        def __req__(self, g):
            return self.__class__(g).ex == self.ex

        def __ne__(self, g):
            return self.ex != self.__class__(g).ex

        def __rne__(self, g):
            return self.__class__(g).ex != self.ex

        def __nonzero__(self):
            return self.ex != 0

    dtype = Expression

    zero  = Expression(0)
    one   = Expression(1)

    rep   = 'EX'

    has_assoc_Ring         = False
    has_assoc_Field        = True

    has_CharacteristicZero = True

    def __init__(self):
        pass

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return a.as_basic()

    def from_sympy(self, a):
        """Convert SymPy's expression to `dtype`. """
        return self.dtype(a)

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_ZZ_sympy(K1, a, K0):
        """Convert a SymPy `Integer` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_QQ_sympy(K1, a, K0):
        """Convert a SymPy `Rational` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_RR_sympy(K1, a, K0):
        """Convert a SymPy `Real` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_PolynomialRing(K1, a, K0):
        """Convert a `DMP` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_FractionField(K1, a, K0):
        """Convert a `DMF` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_ExpressionDomain(K1, a, K0):
        """Convert a `EX` object to `dtype`. """
        return a

    def get_ring(self):
        """Returns a ring associated with `self`. """
        raise DomainError('there is no ring associated with %s' % self)

    def get_field(self):
        """Returns a field associated with `self`. """
        return self

    def is_positive(self, a):
        """Returns True if `a` is positive. """
        return a.ex.as_coeff_terms()[0].is_positive

    def is_negative(self, a):
        """Returns True if `a` is negative. """
        return a.ex.as_coeff_terms()[0].is_negative

    def is_nonpositive(self, a):
        """Returns True if `a` is non-positive. """
        return a.ex.as_coeff_terms()[0].is_nonpositive

    def is_nonnegative(self, a):
        """Returns True if `a` is non-negative. """
        return a.ex.as_coeff_terms()[0].is_nonnegative

    def numer(self, a):
        """Returns the numerator of `a`. """
        return a.numer()

    def denom(self, a):
        """Returns the denominator of `a`. """
        return a.denom()

EX = ExpressionDomain()

