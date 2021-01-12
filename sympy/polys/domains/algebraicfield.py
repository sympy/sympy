"""Implementation of :class:`AlgebraicField` class. """


from sympy.polys.domains.characteristiczero import CharacteristicZero
from sympy.polys.domains.field import Field
from sympy.polys.domains.simpledomain import SimpleDomain
from sympy.polys.polyclasses import ANP
from sympy.polys.polyerrors import CoercionFailed, DomainError, NotAlgebraic, IsomorphismFailed
from sympy.utilities import public

@public
class AlgebraicField(Field, CharacteristicZero, SimpleDomain):
    """A class for representing algebraic number fields. """

    dtype = ANP

    is_AlgebraicField = is_Algebraic = True
    is_Numerical = True

    has_assoc_Ring = False
    has_assoc_Field = True

    def __init__(self, dom, *ext):
        if not dom.is_QQ:
            raise DomainError("ground domain must be a rational field")

        from sympy.polys.numberfields import to_number_field
        if len(ext) == 1 and isinstance(ext[0], tuple):
            self.orig_ext = ext[0][1:]
        else:
            self.orig_ext = ext
        self.ext = to_number_field(ext)
        self.mod = self.ext.minpoly.rep
        self.domain = self.dom = dom

        self.ngens = 1
        self.symbols = self.gens = (self.ext,)
        self.unit = self([dom(1), dom(0)])

        self.zero = self.dtype.zero(self.mod.rep, dom)
        self.one = self.dtype.one(self.mod.rep, dom)

    def new(self, element):
        return self.dtype(element, self.mod.rep, self.dom)

    def __str__(self):
        return str(self.dom) + '<' + str(self.ext) + '>'

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.dom, self.ext))

    def __eq__(self, other):
        """Returns ``True`` if two domains are equivalent. """
        return isinstance(other, AlgebraicField) and \
            self.dtype == other.dtype and self.ext == other.ext

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`. """
        return AlgebraicField(self.dom, *((self.ext,) + extension))

    def to_sympy(self, a):
        """Convert ``a`` to a SymPy object. """
        # Precompute a converter to be reused:
        if not hasattr(self, '_converter'):
            self._converter = _make_converter(self)

        return self._converter(a)

    def from_sympy(self, a):
        """Convert SymPy's expression to ``dtype``. """
        try:
            return self([self.dom.from_sympy(a)])
        except CoercionFailed:
            pass

        from sympy.polys.numberfields import to_number_field

        try:
            return self(to_number_field(a, self.ext).native_coeffs())
        except (NotAlgebraic, IsomorphismFailed):
            raise CoercionFailed(
                "%s is not a valid algebraic number in %s" % (a, self))

    def from_ZZ_python(K1, a, K0):
        """Convert a Python ``int`` object to ``dtype``. """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_python(K1, a, K0):
        """Convert a Python ``Fraction`` object to ``dtype``. """
        return K1(K1.dom.convert(a, K0))

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY ``mpz`` object to ``dtype``. """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY ``mpq`` object to ``dtype``. """
        return K1(K1.dom.convert(a, K0))

    def from_RealField(K1, a, K0):
        """Convert a mpmath ``mpf`` object to ``dtype``. """
        return K1(K1.dom.convert(a, K0))

    def get_ring(self):
        """Returns a ring associated with ``self``. """
        raise DomainError('there is no ring associated with %s' % self)

    def is_positive(self, a):
        """Returns True if ``a`` is positive. """
        return self.dom.is_positive(a.LC())

    def is_negative(self, a):
        """Returns True if ``a`` is negative. """
        return self.dom.is_negative(a.LC())

    def is_nonpositive(self, a):
        """Returns True if ``a`` is non-positive. """
        return self.dom.is_nonpositive(a.LC())

    def is_nonnegative(self, a):
        """Returns True if ``a`` is non-negative. """
        return self.dom.is_nonnegative(a.LC())

    def numer(self, a):
        """Returns numerator of ``a``. """
        return a

    def denom(self, a):
        """Returns denominator of ``a``. """
        return self.one

    def from_AlgebraicField(K1, a, K0):
        """Convert AlgebraicField element 'a' to another AlgebraicField """
        return K1.from_sympy(K0.to_sympy(a))

    def from_GaussianIntegerRing(K1, a, K0):
        """Convert a GaussianInteger element 'a' to ``dtype``. """
        return K1.from_sympy(K0.to_sympy(a))

    def from_GaussianRationalField(K1, a, K0):
        """Convert a GaussianRational element 'a' to ``dtype``. """
        return K1.from_sympy(K0.to_sympy(a))


def _make_converter(K):
    """Construct the converter to convert back to Expr"""
    # Precompute the effect of converting to sympy and expanding expressions
    # like (sqrt(2) + sqrt(3))**2. Asking Expr to do the expansion on every
    # conversion from K to Expr is slow. Here we compute the expansions for
    # each power of the generator and collect together the resulting algebraic
    # terms and the rational coefficients into a matrix.
    from sympy import Add, S

    gen = K.ext.as_expr()
    todom = K.dom.from_sympy

    # We'll let Expr compute the expansions. We won't make any presumptions
    # about what this results in except that it is QQ-linear in some terms
    # that we will call algebraics. The final result will be expressed in
    # terms of those.
    powers = [S.One, gen]
    for n in range(2, K.mod.degree()):
        powers.append((gen * powers[-1]).expand())

    # Collect the rational coefficients and algebraic Expr that can
    # map the ANP coefficients into an expanded sympy expression
    terms = [dict(t.as_coeff_Mul()[::-1] for t in Add.make_args(p)) for p in powers]
    algebraics = set().union(*terms)
    matrix = [[todom(t.get(a, S.Zero)) for t in terms] for a in algebraics]

    # Create a function to do the conversion efficiently:

    def converter(a):
        """Convert a to Expr using converter"""
        from sympy import Add, Mul
        ai = a.rep[::-1]
        tosympy = K.dom.to_sympy
        coeffs_dom = [sum(mij*aj for mij, aj in zip(mi, ai)) for mi in matrix]
        coeffs_sympy = [tosympy(c) for c in coeffs_dom]
        res = Add(*(Mul(c, a) for c, a in zip(coeffs_sympy, algebraics)))
        return res

    return converter
