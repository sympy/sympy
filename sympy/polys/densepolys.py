"""Object--oriented interface to dense polynomial representation. """

from sympy.polys.polyclasses import GenericPoly

class DensePoly(GenericPoly):
    """Dense polynomial over an arbitrary domain. """

    __slots__ = ['rep', 'lev', 'dom', '_hash']

    def __init__(self, rep, dom, lev=None):
        if lev is None:
            rep, lev = dmp_validate(rep)

        self.rep = rep
        self.lev = lev
        self.dom = dom

        self._hash = None

    def __repr__(self):
        return "%s(%s, %s, %s)" % (self.__class__.__name__, self.rep, self.dom)

    def __hash__(self):
        _hash = self._hash

        if _hash is None:
            self._hash = _hash = hash((self.__class__.__name__, repr(self.rep), self.dom))

        return _hash

    def __getstate__(self):
        return (self.rep, self.lev, self.dom, self._hash)

    def __getnewargs__(self):
        return (self.rep, self.lev, self.dom, self._hash)

    def unify(f, g):
        """Unify representations of two multivariate polynomials. """
        if not hasattr(g, '__iter__'):
            if f.lev == g.lev and f.dom == g.dom:
                return f.lev, f.dom, f.per, f.rep, g.rep
            else:
                raise UnificationFailed("can't unify %s with %s" % (f, g))
        else:
            lev, dom, reps = f.lev, f.dom, []

            for gg in g:
                if gg.lev == lev and gg.dom == dom:
                    reps.append(gg.rep)
                else:
                    raise UnificationFailed("can't unify %s with %s" % (f, g))

            return lev, dom, f.per, f.rep, reps

    def per(f, rep, dom=None, lower=False):
        """Create a dense polynomial out of the given representation. """
        lev = f.lev

        if lower:
            if not lev:
                return rep
            else:
                lev -= 1

        if dom is None:
            dom = f.dom

        return DensePoly(rep, dom, lev)

    @classmethod
    def zero(cls, lev, ord, dom):
        """Construct a zero--polynomial with appropriate properties. """
        return cls(dmp_zero(lev), dom, lev)

    @classmethod
    def one(cls, lev, ord, dom):
        """Construct a one--polynomial with appropriate properties. """
        return cls(dmp_one(lev, dom), dom, lev)

    @classmethod
    def from_ground(cls, rep, lev, ord, dom):
        """Create dense representation from an element of the ground domain. """
        return cls(dmp_from_ground(rep, lev, dom), dom, lev)

    @classmethod
    def from_dict(cls, rep, lev, ord, dom):
        """Create dense representation from a ``dict`` with native coefficients. """
        return cls(dmp_from_dict(rep, lev, dom), dom, lev)

    @classmethod
    def from_sympy_dict(cls, rep, lev, ord, dom):
        """Create dense representation from a ``dict`` with SymPy's coefficients. """
        return cls(dmp_from_sympy_dict(rep, lev, dom), dom, lev)

    @classmethod
    def from_list(cls, rep, lev, ord, dom):
        """Create dense representation from a ``list`` with native coefficients. """
        return cls(dmp_from_dict(rep, lev, dom), dom, lev)

    @classmethod
    def from_sympy_list(cls, rep, lev, ord, dom):
        """Create dense representation from a ``list`` with SymPy's coefficients. """
        return cls(dmp_from_sympy_dict(rep, lev, dom), dom, lev)

    def to_ground(f):
        """Convert dense representation to an element of the ground domain. """
        return dmp_to_ground(f.rep, f.lev, f.dom)

    def to_dict(f):
        """Convert dense representation to a ``dict`` with native coefficients. """
        return dmp_to_dict(f.rep, f.lev, f.dom)

    def to_sympy_dict(f):
        """Convert dense representation to a ``dict`` with SymPy's coefficients. """
        return dmp_to_sympy_dict(f.rep, f.lev, f.dom)

    def to_list(f):
        """Convert dense representation to a ``list`` with native coefficients. """
        return dmp_to_dict(f.rep, f.lev, f.dom)

    def to_sympy_list(f):
        """Convert dense representation to a ``list`` with SymPy's coefficients. """
        return dmp_to_sympy_dict(f.rep, f.lev, f.dom)

    def set_domain(f, dom):
        """Set the ground domain in $f$ to ``dom``. """
        if f.dom == dom:
            return f
        else:
            return f.per(dmp_set_domain(f.rep, f.lev, f.dom, dom), dom=dom)

    def ground_to_ring(f):
        """Make the ground domain a ring. """
        return f.set_domain(f.dom.get_ring())

    def ground_to_field(f):
        """Make the ground domain a field. """
        return f.set_domain(f.dom.get_field())

    def ground_to_exact(f):
        """Make the ground domain exact. """
        return f.set_domain(f.dom.get_exact())

    def LC(f):
        """Return the leading coefficient of $f$. """
        return dmp_ground_LC(f.rep, f.lev, f.dom)

    def LM(f):
        """Return the leading monomial of $f$. """
        return dmp_ground_LM(f.rep, f.lev, f.dom)

    def LT(f):
        """Return the leading term of $f$. """
        return dmp_ground_LT(f.rep, f.lev, f.dom)

    def TC(f):
        """Return the trailing coefficient of $f$. """
        return dmp_ground_TC(f.rep, f.lev, f.dom)

    def TM(f):
        """Return the trailing monomial of $f$. """
        return dmp_ground_TM(f.rep, f.lev, f.dom)

    def TT(f):
        """Return the trailing coefficient of $f$. """
        return dmp_ground_TT(f.rep, f.lev, f.dom)

    def EC(f):
        """Return the last non--zero coefficient of $f$. """
        return dmp_ground_EC(f.rep, f.lev, f.dom)

    def EM(f):
        """Return the last non--zero monomial of $f$. """
        return dmp_ground_EM(f.rep, f.lev, f.dom)

    def ET(f):
        """Return the last non--zero coefficient of $f$. """
        return dmp_ground_ET(f.rep, f.lev, f.dom)

    def nth(f, *N):
        """Return $n$--th coefficient of $f$. """
        return dmp_ground_nth(f.rep, N, f.lev, f.dom)

    def coeffs(f):
        """Return all non--zero coefficients of $f$. """
        return dmp_coeffs(f.rep, f.lev, f.dom)

    def monoms(f):
        """Return all non--zero monomials of $f$. """
        return dmp_monoms(f.rep, f.lev, f.dom)

    def terms(f):
        """Return all non--zero terms from $f$. """
        return dmp_terms(f.rep, f.lev, f.dom)

    def all_coeffs(f):
        """Return all coefficients of $f$. """
        return dmp_all_coeffs(f.rep, f.lev, f.dom)

    def all_monoms(f):
        """Return all monomials of $f$. """
        return dmp_all_monoms(f.rep, f.lev, f.dom)

    def all_terms(f):
        """Return all terms of $f$. """
        return dmp_all_terms(f.rep, f.lev, f.dom)

    def degree(f, j=0):
        """Return the degree of $f$ in $x_j$. """
        return dmp_degree_in(f.rep, j, f.lev)

    def degree_list(f):
        """Return the list of degrees of $f$. """
        return dmp_degree_list(f.rep, f.lev)

    def total_degree(f):
        """Return the total degree of $f$. """
        return dmp_total_degree(f.rep, f.lev)

    def deflate(f):
        """Reduce degree of $f$ by mapping $x_i^m$ to $y_i$. """
        J, F = dmp_deflate(f.rep, f.lev, f.dom)
        return J, f.per(F)

    def inflate(f, M):
        """Revert :func:`deflate` by mapping $y_i$ to $x_i^m$. """
        return f.per(dmp_inflate(f.rep, M, f.lev, f.dom))

    def terms_gcd(f):
        """Remove GCD of terms from the polynomial $f$. """
        J, F = dmp_terms_gcd(f.rep, f.lev, f.dom)
        return J, f.per(F)

    def add_ground(f, c):
        """Add an element of the ground domain to $f$. """
        return f.per(dmp_add_ground(f.rep, f.dom.convert(c), f.lev, f.dom))

    def sub_ground(f, c):
        """Subtract an element of the ground domain from $f$. """
        return f.per(dmp_sub_ground(f.rep, f.dom.convert(c), f.lev, f.dom))

    def mul_ground(f, c):
        """Multiply $f$ by an element of the ground domain. """
        return f.per(dmp_mul_ground(f.rep, f.dom.convert(c), f.lev, f.dom))

    def quo_ground(f, c):
        """Quotient of $f$ by an element of the ground domain. """
        return f.per(dmp_quo_ground(f.rep, f.dom.convert(c), f.lev, f.dom))

    def exquo_ground(f, c):
        """Exact quotient of $f$ by an element of the ground domain. """
        return f.per(dmp_exquo_ground(f.rep, f.dom.convert(c), f.lev, f.dom))

    def abs(f):
        """Make all coefficients in $f$ positive. """
        return f.per(dmp_abs(f.rep, f.lev, f.dom))

    def neg(f):
        """Negate all coefficients in $f$. """
        return f.per(dmp_neg(f.rep, f.lev, f.dom))

    def add(f, g):
        """Add two multivariate polynomials $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_add(F, G, lev, dom))

    def sub(f, g):
        """Subtract two multivariate polynomials $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_sub(F, G, lev, dom))

    def mul(f, g):
        """Multiply two multivariate polynomials $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_mul(F, G, lev, dom))

    def sqr(f):
        """Square a multivariate polynomial $f$. """
        return f.per(dmp_sqr(f.rep, f.lev, f.dom))

    def pow(f, n):
        """Raise $f$ to a non--negative power $n$. """
        return f.per(dmp_pow(f.rep, n, f.lev, f.dom))

    def pdiv(f, g):
        """Polynomial pseudo--division of $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        q, r = dmp_pdiv(F, G, lev, dom)
        return per(q), per(r)

    def prem(f, g):
        """Polynomial pseudo--remainder of $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_prem(F, G, lev, dom))

    def pquo(f, g):
        """Polynomial pseudo--quotient of $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_pquo(F, G, lev, dom))

    def pexquo(f, g):
        """Polynomial exact pseudo--quotient of $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_pexquo(F, G, lev, dom))

    def div(f, g):
        """Polynomial division with remainder of $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        q, r = dmp_div(F, G, lev, dom)
        return per(q), per(r)

    def rem(f, g):
        """Compute polynomial remainder of $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_rem(F, G, lev, dom))

    def quo(f, g):
        """Compute polynomial quotient of $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_quo(F, G, lev, dom))

    def exquo(f, g):
        """Compute polynomial exact quotient of $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_exquo(F, G, lev, dom))

    def reduced(f, G):
        """Reduce $f$ modulo a set of polynomials $G$. """
        lev, dom, per, f, G = f.unify(G)
        return per(dmp_reduced(f, G, lev, dom))

    def max_norm(f):
        """Returns maximum norm of $f$. """
        return dmp_max_norm(f.rep, f.lev, f.dom)

    def l1_norm(f):
        """Returns l1 norm of $f$. """
        return dmp_l1_norm(f.rep, f.lev, f.dom)

    def clear_denoms(f, convert=False):
        """Clear denominators in $f$, but keep the ground domain. """
        coeff, F = dmp_clear_denoms(f.rep, f.lev, f.dom, convert=convert)
        return coeff, f.per(F)

    def lift(f):
        """Convert algebraic coefficients to rationals. """
        return f.per(dmp_lift(f.rep, f.lev, f.dom), dom=f.dom.dom)

    def half_gcdex(f, g):
        """Half extended Euclidean algorithm. """
        lev, dom, per, F, G = f.unify(g)
        s, h = dmp_half_gcdex(F, G, dom)
        return per(s), per(h)

    def gcdex(f, g):
        """Extended Euclidean algorithm. """
        lev, dom, per, F, G = f.unify(g)
        s, t, h = dmp_gcdex(F, G, lev, dom)
        return per(s), per(t), per(h)

    def invert(f, g):
        """Invert $f$ modulo $g$, if possible. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_invert(F, G, lev, dom))

    def subresultants(f, g):
        """Compute subresultant PRS sequence of $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        R = dmp_subresultants(F, G, lev, dom)
        return map(per, R)

    def resultant(f, g):
        """Compute resultant of $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_resultant(F, G, lev, dom), lower=True)

    def discriminant(f):
        """Compute discriminant of $f$. """
        return f.per(dmp_discriminant(f.rep, f.lev, f.dom), lower=True)

    def cofactors(f, g):
        """Compute GCD of $f$ and $g$ and their cofactors. """
        lev, dom, per, F, G = f.unify(g)
        h, cff, cfg = dmp_cofactors(F, G, lev, dom)
        return per(h), per(cff), per(cfg)

    def gcd(f, g):
        """Compute polynomial GCD of $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_gcd(F, G, lev, dom))

    def lcm(f, g):
        """Compute polynomial LCM of $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_lcm(F, G, lev, dom))

    def trunc(f, p):
        """Reduce $f$ modulo an element of the ground domain. """
        return f.per(dmp_ground_trunc(f.rep, f.dom.convert(p), f.lev, f.dom))

    def monic(f):
        """Divide all coefficients by the leading coefficient of $f$. """
        return f.per(dmp_ground_monic(f.rep, f.lev, f.dom))

    def content(f):
        """Compute GCD of all coefficients of $f$. """
        return dmp_ground_content(f.rep, f.lev, f.dom)

    def primitive(f):
        """Compute content and the primitive form of $f$. """
        cont, F = dmp_ground_primitive(f.rep, f.lev, f.dom)
        return cont, f.per(F)

    def integrate(f, m=1, j=0):
        """Compute $m$--th order indefinite integral of $f$ in $x_j$. """
        return f.per(dmp_integrate_in(f.rep, m, j, f.lev, f.dom))

    def diff(f, m=1, j=0):
        """Compute $m$--th order derivative of $f$ in $x_j$. """
        return f.per(dmp_diff_in(f.rep, m, j, f.lev, f.dom))

    def eval(f, a, j=0):
        """Evaluate $f$ at the given point $a$ in $x_j$. """
        return f.per(dmp_eval_in(f.rep, f.dom.convert(a), j, f.lev, f.dom), lower=True)

    def mirror(f, j=0):
        """Evaluate efficiently composition $f(-x_j)$. """
        return f.per(dmp_mirror_in(f.rep, j, f.lev, f.dom))

    def scale(f, a, j=0):
        """Evaluate efficiently composition $f(a x_j)$. """
        return f.per(dmp_scale_in(f.rep, f.dom.convert(a), j, f.lev, f.dom))

    def taylor(f, a, j=0):
        """Evaluate efficiently Taylor shift $f(x_j + a)$. """
        return f.per(dmp_taylor_in(f.rep, f.dom.convert(a), j, f.lev, f.dom))

    def transform(f, p, q, j=0):
        """Evaluate functional transformation $q^n \cdot f(p/q)$. """
        lev, dom, per, F, (P, Q) = f.unify((p, q))
        return per(dmp_transform_in(F, P, Q, j, lev, dom))

    def compose(f, g):
        """Compute functional composition of $f$ and $g$. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_compose(F, G, lev, dom))

    def decompose(f):
        """Computes functional decomposition of $f$. """
        return map(f.per, dmp_decompose(f.rep, f.lev, f.dom))

    def sturm(f):
        """Computes the Sturm sequence of $f$. """
        return map(f.per, dmp_sturm(f.rep, f.lev, f.dom))

    def sqf_norm(f):
        """Computes square--free norm of $f$. """
        s, g, r = dmp_sqf_norm(f.rep, f.lev, f.dom)
        return s, f.per(g), f.per(r, dom=f.dom.dom)

    def sqf_part(f):
        """Computes square--free part of $f$. """
        return f.per(dmp_sqf_part(f.rep, f.lev, f.dom))

    def sqf_list(f, all=False):
        """Returns a list of square--free factors of $f$. """
        coeff, factors = dmp_sqf_list(f.rep, f.lev, f.dom, all=all)
        return coeff, [ (f.per(g), k) for g, k in factors ]

    def sqf_list_include(f, all=False):
        """Returns a list of square--free factors of $f$. """
        factors = dmp_sqf_list_include(f.rep, f.lev, f.dom, all=all)
        return [ (f.per(g), k) for g, k in factors ]

    def factor_list(f):
        """Returns a list of irreducible factors of $f$. """
        coeff, factors = dmp_factor_list(f.rep, f.lev, f.dom)
        return coeff, [ (f.per(g), k) for g, k in factors ]

    def factor_list_include(f):
        """Returns a list of irreducible factors of $f$. """
        factors = dmp_factor_list_include(f.rep, f.lev, f.dom)
        return [ (f.per(g), k) for g, k in factors ]

    def real_intervals(f, eps=None, inf=None, sup=None, fast=False, sqf=False):
        """Compute isolating intervals for real roots of $f$. """
        return dmp_real_intervals(f.rep, f.lev, f.dom, eps=eps, inf=inf, sup=sup, fast=fast)

    def complex_intervals(f, eps=None, inf=None, sup=None, fast=False, sqf=False):
        """Compute isolating rectangles for complex roots of $f$. """
        return dmp_complex_intervals(f.rep, f.lev, f.dom, eps=eps, inf=inf, sup=sup, fast=fast)

    def refine_real_root(f, s, t, eps=None, steps=None, fast=False):
        """Refine a real root isolating interval to the given precision. """
        return dmp_refine_real_root(f.rep, s, t, f.lev, f.dom, eps=eps, steps=steps, fast=fast)

    def refine_complex_root(f, s, t, eps=None, steps=None, fast=False):
        """Refine a complex root isolating rectangle to the given precision. """
        return dmp_refine_complex_root(f.rep, s, t, f.lev, f.dom, eps=eps, steps=steps, fast=fast)

    def count_real_roots(f, inf=None, sup=None):
        """Return the number of real roots of $f$ in $[inf, sup]$ interval. """
        return dmp_count_real_roots(f.rep, f.lev, f.dom, inf=inf, sup=sup)

    def count_complex_roots(f, inf=None, sup=None):
        """Return the number of complex roots of $f$ in $[inf, sup]$ rectangle. """
        return dmp_count_complex_roots(f.rep, f.lev, f.dom, inf=inf, sup=sup)

    @property
    def is_zero(f):
        """Returns ``True`` if $f$ is equivalent to zero. """
        return dmp_zero_p(f.rep, f.lev)

    @property
    def is_one(f):
        """Return ``True`` if $f$ is equivalent to one. """
        return dmp_one_p(f.rep, f.lev, f.dom)

    @property
    def is_ground(f):
        """Return ``True`` if $f$ is an element of the ground domain. """
        return dmp_ground_p(f.rep, f.lev)

    @property
    def is_sqf(f):
        """Return ``True`` if $f$ is a square--free polynomial. """
        return dmp_sqf_p(f.rep, f.lev, f.dom)

    @property
    def is_monic(f):
        """Return ``True`` if the leading coefficient of $f$ is one. """
        return dmp_monic_p(f.rep, f.lev, f.dom)

    @property
    def is_primitive(f):
        """Return ``True`` if GCD of coefficients of $f$ is one. """
        return dmp_primitive_p(f.rep, f.lev, f.dom)

    @property
    def is_linear(f):
        """Return ``True`` if $f$ is linear in all its variables. """
        return dmp_linear_p(f.rep, f.lev, f.dom)

    @property
    def is_homogeneous(f):
        """Return ``True`` if $f$ has zero trailing coefficient. """
        return dmp_homogeneous_p(f.rep, f.lev, f.dom)

    def __abs__(f):
        return f.abs()

    def __neg__(f):
        return f.neg()

    def __add__(f, g):
        if not isinstance(g, DensePoly):
            return f.add_ground(g)
        else:
            return f.add(g)

    def __radd__(f, g):
        return f.__add__(g)

    def __sub__(f, g):
        if not isinstance(g, DensePoly):
            return f.sub_ground(g)
        else:
            return f.sub(g)

    def __rsub__(f, g):
        return (-f).__add__(g)

    def __mul__(f, g):
        if not isinstance(g, DensePoly):
            return f.mul_ground(g)
        else:
            return f.mul(g)

    def __rmul__(f, g):
        return f.__mul__(g)

    def __pow__(f, n):
        return f.pow(n)

    def __divmod__(f, g):
        return f.div(g)

    def __mod__(f, g):
        return f.rem(g)

    def __floordiv__(f, g):
        if not isinstance(g, DensePoly):
            return f.exquo_ground(g)
        else:
            return f.exquo(g)

    def __eq__(f, g):
        return isinstance(g, DensePoly) and f.rep == g.rep

    def __ne__(f, g):
        return not f.__eq__(g)

    def __nonzero__(f):
        return not f.is_zero
