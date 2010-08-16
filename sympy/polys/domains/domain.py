"""Implementation of :class:`Domain` class. """

from sympy.core import Basic, sympify

from sympy.polys.polyerrors import (
    UnificationFailed,
    CoercionFailed,
    DomainError,
)

from sympy.polys.domains.groundtypes import python_factorial

class Domain(object):
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

    is_Simple    = False
    is_Composite = False

    has_CharacteristicZero = False

    is_EX = False

    rep   = None
    alias = None

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
        return self.dtype(*args)

    def convert(K1, a, K0=None):
        """Convert an object `a` from `K0` to `K1`. """
        if K0 is not None:
            if K0.alias is not None:
                method = "from_" + K0.alias
            else:
                method = "from_" + K0.__class__.__name__

            _convert = getattr(K1, method)

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

                if K1.is_Numerical and getattr(a, 'is_ground', False):
                    r = a.rep
                    for i in xrange(a.lev + 1):
                        r = r[0]
                    return K1.convert(r)

                a = sympify(a)

                if isinstance(a, Basic):
                    return K1.from_sympy(a)
            except (TypeError, ValueError):
                pass

            raise CoercionFailed("can't convert %s to type %s" % (a, K1))

    def of_type(self, a):
        """Check if `a` is of type `dtype`. """
        return type(a) is type(self.one)

    def __contains__(self, a):
        """Check if `a` belongs to this domain. """
        try:
            self.convert(a)
        except CoercionFailed:
            return False

        return True

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        raise NotImplementedError

    def from_sympy(self, a):
        """Convert a SymPy object to `dtype`. """
        raise NotImplementedError

    def from_FF_python(K1, a, K0):
        """Convert `ModularInteger(int)` to `dtype`. """
        return None

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return None

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return None

    def from_FF_sympy(K1, a, K0):
        """Convert `ModularInteger(Integer)` to `dtype`. """
        return None

    def from_ZZ_sympy(K1, a, K0):
        """Convert a SymPy `Integer` object to `dtype`. """
        return None

    def from_QQ_sympy(K1, a, K0):
        """Convert a SymPy `Rational` object to `dtype`. """
        return None

    def from_FF_gmpy(K1, a, K0):
        """Convert `ModularInteger(mpz)` to `dtype`. """
        return None

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return None

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return None

    def from_RR_sympy(K1, a, K0):
        """Convert a SymPy `Real` object to `dtype`. """
        return None

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
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
        """Returns a maximal domain containg `K0` and `K1`. """
        if gens is not None:
            if (K0.is_Composite and (set(K0.gens) & set(gens))) or (K1.is_Composite and (set(K1.gens) & set(gens))):
                raise UnificationFailed("can't unify %s with %s, given %s generators" % (K0, K1, tuple(gens)))

        if K0 == K1:
            return K0

        if not K0.has_CharacteristicZero:
            if not K1.has_CharacteristicZero:
                if K0.mod == K1.mod and K0.dom == K1.dom:
                    return K0
            elif K1.is_ZZ:
                return K0

            raise UnificationFailed("can't unify %s with %s" % (K0, K1))

        if not K1.has_CharacteristicZero:
            if K0.is_ZZ:
                return K1
            else:
                raise UnificationFailed("can't unify %s with %s" % (K0, K1))

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

        from sympy.polys.domains import EX
        return EX

    def __eq__(self, other):
        """Returns `True` if two domains are equivalent. """
        return self.dtype == other.dtype

    def __ne__(self, other):
        """Returns `False` if two domains are equivalent. """
        return self.dtype != other.dtype

    def map(self, seq):
        """Rersively apply ``self`` to all elements of ``seq``. """
        result = []

        for elt in seq:
            if isinstance(elt, list):
                result.append(self.map(elt))
            else:
                result.append(self(elt))

        return result

    def get_ring(self):
        """Returns a ring associated with `self`. """
        raise DomainError('there is no ring associated with %s' % self)

    def get_field(self):
        """Returns a field associated with `self`. """
        raise DomainError('there is no field associated with %s' % self)

    def get_exact(self):
        """Returns an exact domain associated with `self`. """
        return self

    def float_domain(self):
        return FF

    def complex_domain(self):
        return CC

    def __getitem__(self, gens):
        """The mathematical way do make a polynomial ring. """
        if hasattr(gens, '__iter__'):
            return self.poly_ring(*gens)
        else:
            return self.poly_ring(gens)

    def poly_ring(self, *gens):
        """Returns a polynomial ring, i.e. `K[X]`. """
        from sympy.polys.domains import PolynomialRing
        return PolynomialRing(self, *gens)

    def frac_field(self, *gens):
        """Returns a fraction field, i.e. `K(X)`. """
        from sympy.polys.domains import FractionField
        return FractionField(self, *gens)

    def algebraic_field(self, *extension):
        """Returns an algebraic field, i.e. `K(alpha, ...)`. """
        raise DomainError("can't create algebraic field over %s" % self)

    def inject(self, *gens):
        """Inject generators into this domain. """
        raise NotImplementedError

    def is_zero(self, a):
        """Returns True if `a` is zero. """
        return not a

    def is_one(self, a):
        """Returns True if `a` is one. """
        return a == self.one

    def is_positive(self, a):
        """Returns True if `a` is positive. """
        return a > 0

    def is_negative(self, a):
        """Returns True if `a` is negative. """
        return a < 0

    def is_nonpositive(self, a):
        """Returns True if `a` is non-positive. """
        return a <= 0

    def is_nonnegative(self, a):
        """Returns True if `a` is non-negative. """
        return a >= 0

    def abs(self, a):
        """Absolute value of `a`, implies `__abs__`. """
        return abs(a)

    def neg(self, a):
        """Returns `a` negated, implies `__neg__`. """
        return -a

    def pos(self, a):
        """Returns `a` positive, implies `__pos__`. """
        return +a

    def add(self, a, b):
        """Sum of `a` and `b`, implies `__add__`.  """
        return a + b

    def sub(self, a, b):
        """Difference of `a` and `b`, implies `__sub__`.  """
        return a - b

    def mul(self, a, b):
        """Product of `a` and `b`, implies `__mul__`.  """
        return a * b

    def pow(self, a, b):
        """Raise `a` to power `b`, implies `__pow__`.  """
        return a ** b

    def exquo(self, a, b):
        """Exact quotient of `a` and `b`, implies something. """
        raise NotImplementedError

    def quo(self, a, b):
        """Quotient of `a` and `b`, implies something.  """
        raise NotImplementedError

    def rem(self, a, b):
        """Remainder of `a` and `b`, implies `__mod__`.  """
        raise NotImplementedError

    def div(self, a, b):
        """Division of `a` and `b`, implies something. """
        raise NotImplementedError

    def invert(self, a, b):
        """Returns inversion of `a mod b`, implies something. """
        raise NotImplementedError

    def revert(self, a):
        """Returns `a**(-1)` if possible. """
        raise NotImplementedError

    def numer(self, a):
        """Returns numerator of `a`. """
        raise NotImplementedError

    def denom(self, a):
        """Returns denominator of `a`. """
        raise NotImplementedError

    def gcdex(self, a, b):
        """Extended GCD of `a` and `b`. """
        raise NotImplementedError

    def gcd(self, a, b):
        """Returns GCD of `a` and `b`. """
        raise NotImplementedError

    def lcm(self, a, b):
        """Returns LCM of `a` and `b`. """
        raise NotImplementedError

    def log(self, a, b):
        """Returns b-base logarithm of `a`. """
        raise NotImplementedError

    def sqrt(self, a):
        """Returns square root of `a`. """
        raise NotImplementedError

    def evalf(self, a, **args):
        """Returns numerical approximation of `a`. """
        return self.to_sympy(a).evalf(**args)

    def real(self, a):
        return a

    def imag(self, a):
        return self.zero

    def characteristic(self):
        """Return the characteristic of this domain. """
        raise NotImplementedError('characteristic()')

