"""Object-oriented interface to sparse polynomial representation. """

from sympy.polys.polyclasses import GenericPoly

class SparsePoly(GenericPoly):
    """Sparse polynomial over an arbitrary domain. """

    __slots__ = ['rep', 'lev', 'ord', 'dom', '_hash']

    def __init__(self, rep, ord, dom, lev=None):
        if lev is None:
            rep, lev = smp_validate(rep)

        self.rep = rep
        self.lev = lev
        self.ord = ord
        self.dom = dom

        self._hash = None

    def __repr__(self):
        return "%s(%s, %s, %s)" % (self.__class__.__name__, self.rep, self.ord, self.dom)

    def __hash__(self):
        _hash = self._hash

        if _hash is None:
            self._hash = _hash = hash((self.__class__.__name__, repr(self.rep), self.ord, self.dom))

        return _hash

    def __getstate__(self):
        return (self.rep, self.lev, self.ord, self.dom, self._hash)

    def __getnewargs__(self):
        return (self.rep, self.lev, self.ord, self.dom, self._hash)

    def unify(f, g):
        """Unify representations of two sparse polynomials. """
        if not hasattr(g, '__iter__'):
            if f.lev == g.lev and f.ord == g.ord and f.dom == g.dom:
                return f.lev, f.ord, f.dom, f.per, f.rep, g.rep
            else:
                raise UnificationFailed("can't unify %s with %s" % (f, g))
        else:
            lev, ord, dom, reps = f.lev, f.ord, f.dom, []

            for gg in g:
                if gg.lev == lev and gg.ord == ord and gg.dom == dom:
                    reps.append(gg.rep)
                else:
                    raise UnificationFailed("can't unify %s with %s" % (f, g))

            return lev, ord, dom, f.per, f.rep, reps

    def per(f, rep, ord=None, dom=None, lower=False):
        """Create a sparse polynomial out of the given representation. """
        lev = f.lev

        if lower:
            if not lev:
                return rep
            else:
                lev -= 1

        if dom is None:
            dom = f.dom

        if ord is None:
            ord = f.ord

        return SparsePoly(rep, dom, ord, lev)

    @classmethod
    def zero(cls, lev, ord, dom):
        """Construct a zero-polynomial with appropriate properties. """
        return cls(smp_zero(lev), ord, dom, lev)

    @classmethod
    def one(cls, lev, ord, dom):
        """Construct a one-polynomial with appropriate properties. """
        return cls(smp_one(lev, dom), ord, dom, lev)

    @classmethod
    def from_ground(cls, rep, lev, ord, dom):
        """Create sparse representation from an element of the ground domain. """
        return cls(smp_from_ground(rep, lev, ord, dom), ord, dom, lev)

    @classmethod
    def from_dict(cls, rep, lev, ord, dom):
        """Create sparse representation from a ``dict`` with native coefficients. """
        return cls(smp_from_dict(rep, lev, ord, dom), ord, dom, lev)

    @classmethod
    def from_sympy_dict(cls, rep, lev, ord, dom):
        """Create sparse representation from a ``dict`` with SymPy's coefficients. """
        return cls(smp_from_sympy_dict(rep, lev, ord, dom), ord, dom, lev)

    @classmethod
    def from_list(cls, rep, lev, ord, dom):
        """Create sparse representation from a ``list`` with native coefficients. """
        return cls(smp_from_dict(rep, lev, ord, dom), ord, dom, lev)

    @classmethod
    def from_sympy_list(cls, rep, lev, ord, dom):
        """Create sparse representation from a ``list`` with SymPy's coefficients. """
        return cls(smp_from_sympy_dict(rep, lev, ord, dom), ord, dom, lev)

    def to_ground(f):
        """Convert sparse representation to an element of the ground domain. """
        return smp_to_ground(f.rep, f.lev, f.ord, f.dom)

    def to_dict(f):
        """Convert sparse representation to a ``dict`` with native coefficients. """
        return smp_to_dict(f.rep, f.lev, f.ord, f.dom)

    def to_sympy_dict(f):
        """Convert sparse representation to a ``dict`` with SymPy's coefficients. """
        return smp_to_sympy_dict(f.rep, f.lev, f.ord, f.dom)

    def to_list(f):
        """Convert sparse representation to a ``list`` with native coefficients. """
        return smp_to_dict(f.rep, f.lev, f.ord, f.dom)

    def to_sympy_list(f):
        """Convert sparse representation to a ``list`` with SymPy's coefficients. """
        return smp_to_sympy_dict(f.rep, f.lev, f.ord, f.dom)

    def set_order(f, ord):
        """Set the ordering of monomials in $f$ to ``ord``. """
        if f.ord == ord:
            return f
        else:
            return f.per(smp_set_order(f.rep, f.lev, ord, f.dom), ord=ord)

    def set_domain(f, dom):
        """Set the ground domain in $f$ to ``dom``. """
        if f.dom == dom:
            return f
        else:
            return f.per(smp_set_domain(f.rep, f.lev, f.ord, f.dom, dom), dom=dom)

    def LC(f):
        """Return the leading coefficient of $f$. """
        return smp_ground_LC(f.rep, f.lev, f.ord, f.dom)

    def LM(f):
        """Return the leading monomial of $f$. """
        return smp_ground_LM(f.rep, f.lev, f.ord, f.dom)

    def LT(f):
        """Return the leading term of $f$. """
        return smp_ground_LT(f.rep, f.lev, f.ord, f.dom)

    def TC(f):
        """Return the trailing coefficient of $f$. """
        return smp_ground_TC(f.rep, f.lev, f.ord, f.dom)

    def TM(f):
        """Return the trailing monomial of $f$. """
        return smp_ground_TM(f.rep, f.lev, f.ord, f.dom)

    def TT(f):
        """Return the trailing coefficient of $f$. """
        return smp_ground_TT(f.rep, f.lev, f.ord, f.dom)

    def EC(f):
        """Return the last non-zero coefficient of $f$. """
        return smp_ground_EC(f.rep, f.lev, f.ord, f.dom)

    def EM(f):
        """Return the last non-zero monomial of $f$. """
        return smp_ground_EM(f.rep, f.lev, f.ord, f.dom)

    def ET(f):
        """Return the last non-zero coefficient of $f$. """
        return smp_ground_ET(f.rep, f.lev, f.ord, f.dom)

    def nth(f, *N):
        """Return $n$-th coefficient of $f$. """
        return smp_ground_nth(f.rep, N, f.lev, f.dom)

    def coeffs(f):
        """Return all non-zero coefficients of $f$. """
        return smp_coeffs(f.rep, f.lev, f.ord, f.dom)

    def monoms(f):
        """Return all non-zero monomials of $f$. """
        return smp_monoms(f.rep, f.lev, f.ord, f.dom)

    def terms(f):
        """Return all non-zero terms from $f$. """
        return smp_terms(f.rep, f.lev, f.ord, f.dom)

    def all_coeffs(f):
        """Return all coefficients of $f$. """
        return smp_all_coeffs(f.rep, f.lev, f.ord, f.dom)

    def all_monoms(f):
        """Return all monomials of $f$. """
        return smp_all_monoms(f.rep, f.lev, f.ord, f.dom)

    def all_terms(f):
        """Return all terms of $f$. """
        return smp_all_terms(f.rep, f.lev, f.ord, f.dom)

    def degree(f, j=0):
        """Return the degree of $f$ in $x_j$. """
        return smp_degree(f.rep, j, f.lev)

    def degrees(f):
        """Return the list of degrees of $f$. """
        return smp_degrees(f.rep, f.lev)

    def total_degree(f):
        """Return the total degree of $f$. """
        return smp_total_degree(f.rep, f.lev)

    def deflate(f):
        """Reduce degree of $f$ by mapping $x_i^m$ to $y_i$. """
        M, F = smp_deflate(f.rep, f.lev, f.ord, f.dom)
        return M, f.per(F)

    def inflate(f, M):
        """Revert :func:`deflate` by mapping $y_i$ to $x_i^m$. """
        return f.per(smp_inflate(f.rep, M, f.lev, f.ord, f.dom))

    def terms_gcd(f):
        """Remove GCD of terms from the polynomial $f$. """
        J, F = smp_terms_gcd(f.rep, f.lev, f.ord, f.dom)
        return J, f.per(F)

    def add_ground(f, c):
        """Add an element of the ground domain to $f$. """
        return f.per(smp_add_ground(f.rep, f.dom.convert(c), f.lev, f.ord, f.dom))

    def sub_ground(f, c):
        """Subtract an element of the ground domain from $f$. """
        return f.per(smp_sub_ground(f.rep, f.dom.convert(c), f.lev, f.ord, f.dom))

    def mul_ground(f, c):
        """Multiply $f$ by an element of the ground domain. """
        return f.per(smp_mul_ground(f.rep, f.dom.convert(c), f.lev, f.ord, f.dom))

    def quo_ground(f, c):
        """Quotient of $f$ by an element of the ground domain. """
        return f.per(smp_quo_ground(f.rep, f.dom.convert(c), f.lev, f.ord, f.dom))

    def exquo_ground(f, c):
        """Exact quotient of $f$ by an element of the ground domain. """
        return f.per(smp_exquo_ground(f.rep, f.dom.convert(c), f.lev, f.ord, f.dom))

    def abs(f):
        """Make all coefficients in $f$ positive. """
        return f.per(smp_abs(f.rep, f.lev, f.ord, f.dom))

    def neg(f):
        """Negate all coefficients in $f$. """
        return f.per(smp_neg(f.rep, f.lev, f.ord, f.dom))

    def add(f, g):
        """Add two multivariate polynomials $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_add(F, G, lev, ord, dom))

    def sub(f, g):
        """Subtract two multivariate polynomials $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_sub(F, G, lev, ord, dom))

    def mul(f, g):
        """Multiply two multivariate polynomials $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_mul(F, G, lev, ord, dom))

    def sqr(f):
        """Square a multivariate polynomial $f$. """
        return f.per(smp_sqr(f.rep, f.lev, f.ord, f.dom))

    def pow(f, n):
        """Raise $f$ to a non-negative power $n$. """
        return f.per(smp_pow(f.rep, n, f.lev, f.ord, f.dom))

    def pdiv(f, g):
        """Polynomial pseudo-division of $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        q, r = smp_pdiv(F, G, lev, ord, dom)
        return per(q), per(r)

    def prem(f, g):
        """Polynomial pseudo-remainder of $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_prem(F, G, lev, dom))

    def pquo(f, g):
        """Polynomial pseudo-quotient of $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_pquo(F, G, lev, ord, dom))

    def pexquo(f, g):
        """Polynomial exact pseudo-quotient of $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_pexquo(F, G, lev, ord, dom))

    def div(f, g):
        """Polynomial division with remainder of $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        q, r = smp_div(F, G, lev, ord, dom)
        return per(q), per(r)

    def rem(f, g):
        """Compute polynomial remainder of $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_rem(F, G, lev, ord, dom))

    def quo(f, g):
        """Compute polynomial quotient of $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_quo(F, G, lev, ord, dom))

    def exquo(f, g):
        """Compute polynomial exact quotient of $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_exquo(F, G, lev, ord, dom))

    def reduced(f, G):
        """Reduce $f$ modulo a set of polynomials $G$. """
        lev, ord, dom, per, f, G = f.unify(G)
        return per(smp_reduced(f, G, lev, ord, dom))

    def max_norm(f):
        """Returns maximum norm of $f$. """
        return smp_max_norm(f.rep, f.lev, f.ord, f.dom)

    def l1_norm(f):
        """Returns l1 norm of $f$. """
        return smp_l1_norm(f.rep, f.lev, f.ord, f.dom)

    def clear_denoms(f, convert=False):
        """Clear denominators in $f$, but keep the ground domain. """
        coeff, F = smp_clear_denoms(f.rep, f.lev, f.ord, f.dom, convert=convert)
        return coeff, f.per(F)

    def lift(f):
        """Convert algebraic coefficients to rationals. """
        return f.per(smp_lift(f.rep, f.lev, f.ord, f.dom), dom=f.dom.dom)

    def half_gcdex(f, g):
        """Half extended Euclidean algorithm. """
        lev, ord, dom, per, F, G = f.unify(g)
        s, h = smp_half_gcdex(F, G, ord, dom)
        return per(s), per(h)

    def gcdex(f, g):
        """Extended Euclidean algorithm. """
        lev, ord, dom, per, F, G = f.unify(g)
        s, t, h = smp_gcdex(F, G, lev, ord, dom)
        return per(s), per(t), per(h)

    def invert(f, g):
        """Invert $f$ modulo $g$, if possible. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_invert(F, G, lev, ord, dom))

    def subresultants(f, g):
        """Compute subresultant PRS sequence of $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        R = smp_subresultants(F, G, lev, ord, dom)
        return map(per, R)

    def resultant(f, g):
        """Compute resultant of $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_resultant(F, G, lev, ord, dom), lower=True)

    def discriminant(f):
        """Compute discriminant of $f$. """
        return f.per(smp_discriminant(f.rep, f.lev, f.ord, f.dom), lower=True)

    def cofactors(f, g):
        """Compute GCD of $f$ and $g$ and their cofactors. """
        lev, ord, dom, per, F, G = f.unify(g)
        h, cff, cfg = smp_cofactors(F, G, lev, ord, dom)
        return per(h), per(cff), per(cfg)

    def gcd(f, g):
        """Compute polynomial GCD of $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_gcd(F, G, lev, ord, dom))

    def lcm(f, g):
        """Compute polynomial LCM of $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_lcm(F, G, lev, ord, dom))

    def trunc(f, p):
        """Reduce $f$ modulo an element of the ground domain. """
        return f.per(smp_ground_trunc(f.rep, f.dom.convert(p), f.lev, f.ord, f.dom))

    def monic(f):
        """Divide all coefficients by the leading coefficient of $f$. """
        return f.per(smp_ground_monic(f.rep, f.lev, f.ord, f.dom))

    def content(f):
        """Compute GCD of all coefficients of $f$. """
        return smp_ground_content(f.rep, f.lev, f.ord, f.dom)

    def primitive(f):
        """Compute content and the primitive form of $f$. """
        cont, F = smp_ground_primitive(f.rep, f.lev, f.ord, f.dom)
        return cont, f.per(F)

    def integrate(f, m=1, j=0):
        """Compute $m$-th order indefinite integral of $f$ in $x_j$. """
        return f.per(smp_integrate_in(f.rep, m, j, f.lev, f.ord, f.dom))

    def diff(f, m=1, j=0):
        """Compute $m$-th order derivative of $f$ in $x_j$. """
        return f.per(smp_diff_in(f.rep, m, j, f.lev, f.ord, f.dom))

    def eval(f, a, j=0):
        """Evaluate $f$ at the given point $a$ in $x_j$. """
        return f.per(smp_eval_in(f.rep, f.dom.convert(a), j, f.lev, f.ord, f.dom), lower=True)

    def mirror(f, j=0):
        """Evaluate efficiently composition $f(-x_j)$. """
        return f.per(smp_mirror_in(f.rep, j, f.lev, f.ord, f.dom))

    def scale(f, a, j=0):
        """Evaluate efficiently composition $f(a x_j)$. """
        return f.per(smp_scale_in(f.rep, f.dom.convert(a), j, f.lev, f.ord, f.dom))

    def taylor(f, a, j=0):
        """Evaluate efficiently Taylor shift $f(x_j + a)$. """
        return f.per(smp_taylor_in(f.rep, f.dom.convert(a), j, f.lev, f.ord, f.dom))

    def transform(f, p, q, j=0):
        """Evaluate functional transformation $q^n \cdot f(p/q)$. """
        lev, ord, dom, per, F, (P, Q) = f.unify((p, q))
        return per(smp_transform_in(F, P, Q, j, lev, ord, dom))

    def compose(f, g):
        """Compute functional composition of $f$ and $g$. """
        lev, ord, dom, per, F, G = f.unify(g)
        return per(smp_compose(F, G, lev, ord, dom))

    def decompose(f):
        """Computes functional decomposition of $f$. """
        return map(f.per, smp_decompose(f.rep, f.lev, f.ord, f.dom))

    def sturm(f):
        """Computes the Sturm sequence of $f$. """
        return map(f.per, smp_sturm(f.rep, f.lev, f.ord, f.dom))

    def sqf_norm(f):
        """Compute square-free norm of $f$. """
        s, g, r = smp_sqf_norm(f.rep, f.lev, f.ord, f.dom)
        return s, f.per(g), f.per(r, dom=f.dom.dom)

    def sqf_part(f):
        """Compute square-free part of $f$. """
        return f.per(smp_sqf_part(f.rep, f.lev, f.ord, f.dom))

    def sqf_list(f, all=False, include=False):
        """Return a list of square-free factors of $f$. """
        result = smp_sqf_list(f.rep, f.lev, f.ord, f.dom, all=all, include=include)
        return f._perify_factors(result, include)

    def factor_list(f, include=False):
        """Return a list of irreducible factors of $f$. """
        result = smp_factor_list(f.rep, f.lev, f.ord, f.dom, include=include)
        return f._perify_factors(f.per, result, include)

    def real_intervals(f, eps=None, inf=None, sup=None, fast=False, sqf=False):
        """Compute isolating intervals for real roots of $f$. """
        return smp_real_intervals(f.rep, f.lev, f.ord, f.dom, eps=eps, inf=inf, sup=sup, fast=fast, sqf=sqf)

    def complex_intervals(f, eps=None, inf=None, sup=None, fast=False, sqf=False):
        """Compute isolating rectangles for complex roots of $f$. """
        return smp_complex_intervals(f.rep, f.lev, f.ord, f.dom, eps=eps, inf=inf, sup=sup, fast=fast, sqf=sqf)

    def refine_real_root(f, s, t, eps=None, steps=None, fast=False, sqf=False):
        """Refine a real root isolating interval to the given precision. """
        return smp_refine_real_root(f.rep, s, t, f.lev, f.ord, f.dom, eps=eps, steps=steps, fast=fast, sqf=sqf)

    def refine_complex_root(f, s, t, eps=None, steps=None, fast=False, sqf=False):
        """Refine a complex root isolating rectangle to the given precision. """
        return smp_refine_complex_root(f.rep, s, t, f.lev, f.ord, f.dom, eps=eps, steps=steps, fast=fast, sqf=sqf)

    def count_real_roots(f, inf=None, sup=None):
        """Return the number of real roots of $f$ in $[inf, sup]$ interval. """
        return smp_count_real_roots(f.rep, f.lev, f.ord, f.dom, inf=inf, sup=sup)

    def count_complex_roots(f, inf=None, sup=None):
        """Return the number of complex roots of $f$ in $[inf, sup]$ rectangle. """
        return smp_count_complex_roots(f.rep, f.lev, f.ord, f.dom, inf=inf, sup=sup)

    @property
    def is_zero(f):
        """Returns ``True`` if $f$ is equivalent to zero. """
        return smp_zero_p(f.rep, f.lev)

    @property
    def is_one(f):
        """Return ``True`` if $f$ is equivalent to one. """
        return smp_one_p(f.rep, f.lev, f.dom)

    @property
    def is_ground(f):
        """Return ``True`` if $f$ is an element of the ground domain. """
        return smp_ground_p(f.rep, f.lev)

    @property
    def is_sqf(f):
        """Return ``True`` if $f$ is a square-free polynomial. """
        return smp_sqf_p(f.rep, f.lev, f.ord, f.dom)

    @property
    def is_monic(f):
        """Return ``True`` if the leading coefficient of $f$ is one. """
        return smp_monic_p(f.rep, f.lev, f.ord, f.dom)

    @property
    def is_primitive(f):
        """Return ``True`` if GCD of coefficients of $f$ is one. """
        return smp_primitive_p(f.rep, f.lev, f.ord, f.dom)

    @property
    def is_linear(f):
        """Return ``True`` if $f$ is linear in all its variables. """
        return smp_linear_p(f.rep, f.lev, f.ord, f.dom)

    @property
    def is_homogeneous(f):
        """Return ``True`` if $f$ has zero trailing coefficient. """
        return smp_homogeneous_p(f.rep, f.lev, f.ord, f.dom)

    def __abs__(f):
        return f.abs()

    def __neg__(f):
        return f.neg()

    def __add__(f, g):
        if not isinstance(g, SparsePoly):
            return f.add_ground(g)
        else:
            return f.add(g)

    def __radd__(f, g):
        return f.__add__(g)

    def __sub__(f, g):
        if not isinstance(g, SparsePoly):
            return f.sub_ground(g)
        else:
            return f.sub(g)

    def __rsub__(f, g):
        return (-f).__add__(g)

    def __mul__(f, g):
        if not isinstance(g, SparsePoly):
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
        if not isinstance(g, SparsePoly):
            return f.exquo_ground(g)
        else:
            return f.exquo(g)

    def __eq__(f, g):
        return isinstance(g, SparsePoly) and f.rep == g.rep

    def __ne__(f, g):
        return not f.__eq__(g)

    def __nonzero__(f):
        return not f.is_zero
