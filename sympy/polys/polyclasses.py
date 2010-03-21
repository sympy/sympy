"""OO layer for several polynomial representations. """

from sympy.utilities import any, all

from sympy.polys.densebasic import (
    dmp_validate,
    dup_normal, dmp_normal,
    dup_convert, dmp_convert,
    dup_from_sympy, dmp_from_sympy,
    dup_strip, dmp_strip,
    dup_degree, dmp_degree_in,
    dmp_degree_list,
    dmp_negative_p, dmp_positive_p,
    dup_LC, dmp_ground_LC,
    dup_TC, dmp_ground_TC,
    dup_nth, dmp_ground_nth,
    dmp_zero, dmp_one, dmp_ground,
    dmp_zero_p, dmp_one_p, dmp_ground_p,
    dup_from_dict, dmp_from_dict,
    dup_to_raw_dict, dmp_to_dict,
    dup_deflate, dmp_deflate,
    dup_terms_gcd, dmp_terms_gcd,
    dmp_list_terms,
)

from sympy.polys.densearith import (
    dup_add_term, dmp_add_term,
    dup_sub_term, dmp_sub_term,
    dup_mul_term, dmp_mul_term,
    dup_mul_ground, dmp_mul_ground,
    dup_quo_ground, dmp_quo_ground,
    dup_exquo_ground, dmp_exquo_ground,
    dup_abs, dmp_abs,
    dup_neg, dmp_neg,
    dup_add, dmp_add,
    dup_sub, dmp_sub,
    dup_mul, dmp_mul,
    dup_sqr, dmp_sqr,
    dup_pow, dmp_pow,
    dup_pdiv, dmp_pdiv,
    dup_prem, dmp_prem,
    dup_pquo, dmp_pquo,
    dup_pexquo, dmp_pexquo,
    dup_div, dmp_div,
    dup_rem, dmp_rem,
    dup_quo, dmp_quo,
    dup_exquo, dmp_exquo,
    dmp_add_mul, dmp_sub_mul,
    dup_max_norm, dmp_max_norm,
    dup_l1_norm, dmp_l1_norm,
)

from sympy.polys.densetools import (
    dup_ground_to_ring, dmp_ground_to_ring,
    dup_integrate, dmp_integrate_in,
    dup_diff, dmp_diff_in,
    dup_eval, dmp_eval_in,
    dup_half_gcdex, dup_gcdex, dup_invert,
    dup_subresultants, dmp_subresultants,
    dup_resultant, dmp_resultant,
    dup_discriminant, dmp_discriminant,
    dup_inner_gcd, dmp_inner_gcd,
    dup_gcd, dmp_gcd,
    dup_lcm, dmp_lcm,
    dup_trunc, dmp_ground_trunc,
    dup_content, dmp_ground_content,
    dup_primitive, dmp_ground_primitive,
    dup_monic, dmp_ground_monic,
    dup_sqf_p, dmp_sqf_p,
    dup_sqf_norm, dmp_sqf_norm,
    dup_sqf_part, dmp_sqf_part,
    dup_sqf_list, dup_sqf_list_include,
    dmp_sqf_list, dmp_sqf_list_include,
    dup_compose, dmp_compose,
    dup_decompose,
    dup_sturm,
    dmp_lift,
    dup_isolate_real_roots_sqf,
    dup_isolate_real_roots,
)

from sympy.polys.factortools import (
    dup_factor_list, dup_factor_list_include,
    dmp_factor_list, dmp_factor_list_include,
)

from sympy.polys.galoistools import (
    gf_degree,
    gf_int,
    gf_LC, gf_TC,
    gf_from_dict, gf_to_dict,
    gf_from_int_poly, gf_to_int_poly,
    gf_trunc, gf_normal, gf_convert,
    gf_neg,
    gf_add_ground, gf_sub_ground,
    gf_mul_ground, gf_exquo_ground,
    gf_add, gf_sub, gf_mul, gf_sqr, gf_pow,
    gf_div, gf_rem, gf_quo, gf_exquo,
    gf_gcd, gf_gcdex,
    gf_monic, gf_diff, gf_eval, gf_compose,
    gf_sqf_p, gf_sqf_part, gf_sqf_list,
    gf_factor, gf_irreducible_p,
)

from sympy.polys.groebnertools import (
    sdp_LC, sdp_LM, sdp_LT,
    sdp_coeffs, sdp_monoms,
    sdp_sort, sdp_strip,
    sdp_normal, #sdp_convert,
    sdp_from_dict, sdp_to_dict,
    sdp_one_p, sdp_one,
    sdp_abs, sdp_neg,
    sdp_add_term, sdp_sub_term, sdp_mul_term,
    sdp_add, sdp_sub, sdp_mul, sdp_sqr, sdp_pow,
    sdp_monic, sdp_content, sdp_primitive,
    sdp_div, sdp_rem, sdp_quo, sdp_exquo,
    sdp_lcm, sdp_gcd,
)

from sympy.polys.polyerrors import (
    UnificationFailed,
    PolynomialError,
    DomainError,
)

def init_normal_GFP(rep, mod, dom):
    return GFP(gf_normal(rep, mod, dom), mod, dom)

class GFP(object):
    """Univariate Polynomials over Galois Fields. """

    __slots__ = ['rep', 'mod', 'lev', 'dom', 'sym']

    def __init__(self, rep, mod, dom, symmetric=None):
        if not dom.is_ZZ:
            raise DomainError("only ZZ domains allowed in GFP")

        if type(rep) is dict:
            self.rep = gf_from_dict(rep, mod, dom)
        else:
            if type(rep) is not list:
                self.rep = gf_normal([rep], mod, dom)
            else:
                self.rep = gf_trunc(rep, mod)

        self.mod = mod
        self.lev = 0
        self.dom = dom

        if symmetric is not None:
            self.sym = symmetric
        else:
            self.sym = True

    def __repr__(f):
        return "%s(%s, %s, %s)" % (f.__class__.__name__, f.rep, f.mod, f.dom)

    def __hash__(f):
        return hash((f.__class__.__name__, repr(f.rep), f.mod, f.dom))

    def __getstate__(self):
        return (self.rep, self.mod, self.dom)

    def __getnewargs__(self):
        return (self.rep, self.mod, self.dom)

    def unify(f, g):
        """Unify representations of two GFP polynomials. """
        if not isinstance(g, GFP) or f.mod != g.mod:
            raise UnificationFailed("can't unify %s with %s" % (f, g))

        sym = max(f.sym, g.sym)

        if f.dom == g.dom:
            return f.mod, f.dom, f.per, f.rep, g.rep
        else:
            mod, dom = f.mod, f.dom.unify(g.dom)

            F = gf_convert(f.rep, mod, f.dom, dom)
            G = gf_convert(g.rep, mod, g.dom, dom)

            def per(rep, mod=mod, dom=dom, sym=sym):
                return GFP(rep, mod, dom, sym)

        return mod, dom, per, F, G

    def per(f, rep):
        """Create a GFP out of the given representation. """
        return GFP(rep, f.mod, f.dom, f.sym)

    def to_dict(f):
        """Convert `f` to a dict representation with native coefficients. """
        rep = gf_to_dict(f.rep, f.mod, f.sym)

        for k, v in dict(rep).iteritems():
            rep[(k,)] = v
            del rep[k]

        return rep

    def to_sympy_dict(f):
        """Convert `f` to a dict representation with SymPy coefficients. """
        rep = gf_to_dict(f.rep, f.mod, f.sym)

        for k, v in dict(rep).iteritems():
            rep[(k,)] = f.dom.to_sympy(v)
            del rep[k]

        return rep

    def to_field(f):
        """Make the ground domain a field. """
        return f

    @classmethod
    def zero(cls, dom, mod):
        return GFP(0, mod, dom)

    @classmethod
    def one(cls, dom):
        return GFP(1, mod, dom)

    def trunc(f, mod):
        """Reduce `f` using new modulus. """
        if mod == f.mod:
            return f
        else:
            return GFP(gf_trunc(f.rep, mod), mod, f.dom, f.sym)

    def convert(f, dom):
        """Convert the ground domain of `f`. """
        if f.dom == dom:
            return f
        elif dom.is_ZZ:
            return GFP(gf_convert(f.rep, f.mod, f.dom, dom), f.mod, dom, f.sym)
        else:
            raise DomainError("can't convert GFP ground domain to %s" % dom)

    def coeffs(f):
        """Returns all non-zero coefficients from `f`. """
        if not f:
            return [f.dom.zero]
        elif not f.sym:
            return [ c for c in f.rep if c ]
        else:
            return [ gf_int(c, f.mod) for c in f.rep if c ]

    def monoms(f):
        """Returns all non-zero monomials from `f`. """
        n = gf_degree(f.rep)

        if n < 0:
            return [(0,)]
        else:
            return [ (n-i,) for i, c in enumerate(f.rep) if c ]

    def terms(f):
        """Returns all non-zero terms from `f`. """
        n = gf_degree(f.rep)

        if n < 0:
            return [((0,), f.dom.zero)]
        elif not f.sym:
            return [ ((n-i,), c) for i, c in enumerate(f.rep) if c ]
        else:
            return [ ((n-i,), gf_int(c, f.mod)) for i, c in enumerate(f.rep) if c ]

    def all_coeffs(f):
        """Returns all coefficients from `f`. """
        if not f:
            return [f.dom.zero]
        elif not f.sym:
            return [ c for c in f.rep ]
        else:
            return [ gf_int(c, f.mod) for c in f.rep ]

    def all_monoms(f):
        """Returns all monomials from `f`. """
        if n < 0:
            return [((0,), f.dom.zero)]
        else:
            return [ (n-i,) for i, c in enumerate(f.rep) ]

    def all_terms(f):
        """Returns all terms from a `f`. """
        n = gf_degree(f.rep)

        if n < 0:
            return [((0,), f.dom.zero)]
        elif not f.sym:
            return [ ((n-i,), c) for i, c in enumerate(f.rep) ]
        else:
            return [ ((n-i,), gf_int(c, f.mod)) for i, c in enumerate(f.rep) ]

    def deflate(f):
        """Reduce degree of `f` by mapping `x**m` to `y`. """
        j, F = dup_deflate(f.rep, f.dom)
        return j, f.per(F)

    def terms_gcd(f):
        """Remove GCD of terms from the polynomial `f`. """
        j, F = dup_terms_gcd(f.rep, f.dom)
        return j, f.per(F)

    def add_ground(f, c):
        """Add an element of the ground domain to `f`. """
        return f.per(gf_add_ground(f.rep, f.dom.convert(c), f.mod, f.dom))

    def sub_ground(f, c):
        """Subtract an element of the ground domain from `f`. """
        return f.per(gf_sub_ground(f.rep, f.dom.convert(c), f.mod, f.dom))

    def mul_ground(f, c):
        """Multiply `f` by an element of the ground domain. """
        return f.per(gf_mul_ground(f.rep, f.dom.convert(c), f.mod, f.dom))

    def exquo_ground(f, c):
        """Divide `f` by an element of the ground domain. """
        return f.per(gf_exquo_ground(f.rep, f.dom.convert(c), f.mod, f.dom))

    def neg(f):
        """Negate all cefficients in `f`. """
        return f.per(gf_neg(f.rep, f.dom))

    def add(f, g):
        """Add two univariate polynomials `f` and `g`. """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_add(F, G, mod, dom))

    def sub(f, g):
        """Subtract two univariate polynomials `f` and `g`. """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_sub(F, G, mod, dom))

    def mul(f, g):
        """Multiply two univariate polynomials `f` and `g`. """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_mul(F, G, mod, dom))

    def sqr(f):
        """Square a univariate polynomial `f`. """
        return f.per(gf_sqr(f.rep, f.mod, f.dom))

    def pow(f, n):
        """Raise `f` to a non-negative power `n`. """
        if isinstance(n, int):
            return f.per(gf_pow(f.rep, n, f.mod, f.dom))
        else:
            raise TypeError("`int` expected, got %s" % type(n))

    def div(f, g):
        """Polynomial division with remainder of `f` and `g`. """
        mod, dom, per, F, G = f.unify(g)
        q, r = gf_div(F, G, mod, dom)
        return per(q), per(r)

    def rem(f, g):
        """Computes polynomial remainder in of `f` and `g`. """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_rem(F, G, mod, dom))

    def quo(f, g):
        """Computes polynomial quotient in of `f` and `g`. """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_quo(F, G, mod, dom))

    def exquo(f, g):
        """Computes polynomial exact quotient in of `f` and `g`. """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_exquo(F, G, mod, dom))

    def degree(f):
        """Returns the leading degree of `f`. """
        return gf_degree(f.rep)

    def LC(f):
        """Returns the leading coefficent of `f`. """
        return gf_LC(f.rep, f.dom)

    def TC(f):
        """Returns the trailing coefficent of `f`. """
        return gf_TC(f.rep, f.dom)

    def gcdex(f, g):
        """Extended Euclidean algorithm. """
        mod, dom, per, F, G = f.unify(g)
        s, t, h = gf_gcdex(F, G, mod, dom)
        return per(s), per(t), per(h)

    def gcd(f, g):
        """Returns polynomial GCD of `f` and `g`. """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_gcd(F, G, dom, mod))

    def monic(f):
        """Divides all coefficients by `LC(f)`. """
        return f.per(gf_monic(f.rep, f.mod, f.dom)[1])

    def diff(f, m=1):
        """Computes partial derivative of `f`. """
        if isinstance(m, int):
            if not m:
                return f

            rep = f.rep

            for i in xrange(0, m):
                rep = gf_diff(f.rep, f.mod, f.dom)

            return f.per(rep)
        else:
            raise TypeError("`int` expected, got %s" % type(m))

    def eval(f, a):
        """Evaluates `f` at the given point `a`. """
        return gf_eval(f.rep, f.dom.convert(a), f.mod, f.dom)

    def compose(f, g):
        """Computes functional composition of `f` and `g`. """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_compose(F, G, mod, dom))

    def sqf_part(f):
        """Computes square-free part of `f`. """
        return f.per(gf_sqf_part(f.rep, f.mod, f.dom))

    def sqf_list(f, **args):
        """Returns a list of square-free factors of `f`. """
        coeff, factors = gf_sqf_list(f.rep, f.mod, f.dom, **args)
        return coeff, [ (f.per(g), k) for g, k in factors ]

    def factor_list(f, **args):
        """Returns a list of irreducible factors of `f`. """
        coeff, factors = gf_factor(f.rep, f.mod, f.dom, **args)
        return coeff, [ (f.per(g), k) for g, k in factors ]

    @property
    def is_zero(f):
        """Returns `True` if `f` is a zero polynomial. """
        return not f

    @property
    def is_one(f):
        """Returns `True` if `f` is a unit polynomial. """
        return f.rep == [f.dom.one]

    @property
    def is_ground(f):
        """Returns `True` if `f` is an element of the ground domain. """
        return len(f.rep) <= 1

    @property
    def is_sqf(f):
        """Returns `True` if `f` is a square-free polynomial. """
        return gf_sqf_p(f.rep, f.mod, f.dom)

    @property
    def is_monic(f):
        """Returns `True` if the leading coefficient of `f` is one. """
        return f.dom.is_one(gf_LC(f.rep, f.dom))

    @property
    def is_irreducible(f):
        """Returns `True` if `f` has no factors over its domain. """
        return gf_irreducible_p(f.rep, f.mod, f.dom)

    def __neg__(f):
        return f.neg()

    def __add__(f, g):
        if isinstance(g, GFP):
            return f.add(g)
        else:
            try:
                return f.add_ground(g)
            except TypeError:
                return NotImplemented

    def __radd__(f, g):
        return f.__add__(g)

    def __sub__(f, g):
        if isinstance(g, GFP):
            return f.sub(g)
        else:
            try:
                return f.sub_ground(g)
            except TypeError:
                return NotImplemented

    def __rsub__(f, g):
        return (-f).__add__(g)

    def __mul__(f, g):
        if isinstance(g, GFP):
            return f.mul(g)
        else:
            try:
                return f.mul_ground(g)
            except TypeError:
                return NotImplemented

    def __rmul__(f, g):
        return f.__mul__(g)

    def __pow__(f, n):
        return f.pow(n)

    def __divmod__(f, g):
        return f.div(g)

    def __mod__(f, g):
        return f.rem(g)

    def __floordiv__(f, g):
        if isinstance(g, GFP):
            return f.exquo(g)
        else:
            try:
                return f.exquo_ground(g)
            except TypeError:
                return NotImplemented

    def __eq__(f, g):
        try:
            _, _, _, F, G = f.unify(g)

            return F == G
        except UnificationFailed:
            return False

    def __ne__(f, g):
        try:
            _, _, _, F, G = f.unify(g)

            return F != G
        except UnificationFailed:
            return True

    def __nonzero__(f):
        return bool(f.rep)

def init_normal_DUP(rep, dom):
    return DUP(dup_normal(rep, dom), dom)

class DUP(object):
    """Dense Univariate Polynomials over `K`. """

    __slots__ = ['rep', 'lev', 'dom']

    def __init__(self, rep, dom):
        if type(rep) is dict:
            self.rep = dup_from_dict(rep, dom)
        else:
            if type(rep) is not list:
                rep = [dom.convert(rep)]

            self.rep = dup_strip(rep)

        self.lev = 0
        self.dom = dom

    def __repr__(f):
        return "%s(%s, %s)" % (f.__class__.__name__, f.rep, f.dom)

    def __hash__(f):
        return hash((f.__class__.__name__, repr(f.rep), f.dom))

    def __getstate__(self):
        return (self.rep, self.dom)

    def __getnewargs__(self):
        return (self.rep, self.dom)

    def unify(f, g):
        """Unify representations of two univariate polynomials. """
        if not isinstance(g, DUP):
            raise UnificationFailed("can't unify %s with %s" % (f, g))

        if f.dom == g.dom:
            return f.dom, f.per, f.rep, g.rep
        else:
            dom = f.dom.unify(g.dom)

            F = dup_convert(f.rep, f.dom, dom)
            G = dup_convert(g.rep, g.dom, dom)

            def per(rep, dom=dom):
                return DUP(rep, dom)

        return dom, per, F, G

    def per(f, rep):
        """Create a DUP out of the given representation. """
        return DUP(rep, f.dom)

    @classmethod
    def zero(cls, dom):
        return DUP(0, dom)

    @classmethod
    def one(cls, dom):
        return DUP(1, dom)

    def to_dict(f):
        """Convert `f` to a dict representation with native coefficients. """
        return dup_to_raw_dict(f.rep)

    def to_sympy_dict(f):
        """Convert `f` to a dict representation with SymPy coefficients. """
        rep = dup_to_raw_dict(f.rep)

        for k, v in rep.iteritems():
            rep[k] = f.dom.to_sympy(v)

        return rep

    def to_ring(f):
        """Make the ground domain a ring. """
        return f.convert(f.dom.get_ring())

    def to_field(f):
        """Make the ground domain a field. """
        return f.convert(f.dom.get_field())

    def convert(f, dom):
        """Convert the ground domain of `f`. """
        if f.dom == dom:
            return f
        else:
            return DUP(dup_convert(f.rep, f.dom, dom), dom)

    def coeffs(f):
        """Returns all non-zero coefficients from `f`. """
        if not f:
            return [K.zero]
        else:
            return [ c for c in f.rep if c ]

    def monoms(f):
        """Returns all non-zero monomials from `f`. """
        n = dup_degree(f.rep)

        if n < 0:
            return [(0,)]
        else:
            return [ (n-i,) for i, c in enumerate(f.rep) if c ]

    def terms(f):
        """Returns all non-zero terms from `f`. """
        n = dup_degree(f.rep)

        if n < 0:
            return [((0,), K.zero)]
        else:
            return [ ((n-i,), c) for i, c in enumerate(f.rep) if c ]

    def all_coeffs(f):
        """Returns all coefficients from `f`. """
        if not f:
            return [K.zero]
        else:
            return [ c for c in f.rep ]

    def all_monoms(f):
        """Returns all monomials from `f`. """
        n = dup_degree(f.rep)

        if n < 0:
            return [(0,)]
        else:
            return [ (n-i,) for i, c in enumerate(f.rep) ]

    def all_terms(f):
        """Returns all terms from a `f`. """
        n = dup_degree(f.rep)

        if n < 0:
            return [((0,), K.zero)]
        else:
            return [ ((n-i,), c) for i, c in enumerate(f.rep) ]

    def add_term(f, c, i=0):
        """Add a single term of the form `c*x**i` to `f`. """
        if isinstance(i, int):
            return f.per(dup_add_term(f.rep, f.dom.convert(c), i, f.dom))
        else:
            raise TypeError("`int` expected, got %s" % type(i))

    def sub_term(f, c, i=0):
        """Subtract a single term of the form `c*x**i` from `f`. """
        if isinstance(i, int):
            return f.per(dup_sub_term(f.rep, f.dom.convert(c), i, f.dom))
        else:
            raise TypeError("`int` expected, got %s" % type(i))

    def mul_term(f, c, i=0):
        """Multiply `f` by a single term of the form `c*x**i`. """
        if isinstance(i, int):
            return f.per(dup_mul_term(f.rep, f.dom.convert(c), i, f.dom))
        else:
            raise TypeError("`int` expected, got %s" % type(i))

    def mul_ground(f, c):
        """Multiply `f` by an element of the ground domain. """
        return f.per(dup_mul_ground(f.rep, f.dom.convert(c), f.dom))

    def quo_ground(f, c):
        """Quotient of `f` by an element of the ground domain. """
        return f.per(dup_quo_ground(f.rep, f.dom.convert(c), f.dom))

    def exquo_ground(f, c):
        """Exact quotient of `f` by an element of the ground domain. """
        return f.per(dup_exquo_ground(f.rep, f.dom.convert(c), f.dom))

    def abs(f):
        """Make all coefficients in `f` positive. """
        return f.per(dup_abs(f.rep, f.dom))

    def neg(f):
        """Negate all cefficients in `f`. """
        return f.per(dup_neg(f.rep, f.dom))

    def add(f, g):
        """Add two univariate polynomials `f` and `g`. """
        dom, per, F, G = f.unify(g)
        return per(dup_add(F, G, dom))

    def sub(f, g):
        """Subtract two univariate polynomials `f` and `g`. """
        dom, per, F, G = f.unify(g)
        return per(dup_sub(F, G, dom))

    def mul(f, g):
        """Multiply two univariate polynomials `f` and `g`. """
        dom, per, F, G = f.unify(g)
        return per(dup_mul(F, G, dom))

    def sqr(f):
        """Square a univariate polynomial `f`. """
        return f.per(dup_sqr(f.rep, f.dom))

    def pow(f, n):
        """Raise `f` to a non-negative power `n`. """
        if isinstance(n, int):
            return f.per(dup_pow(f.rep, n, f.dom))
        else:
            raise TypeError("`int` expected, got %s" % type(n))

    def pdiv(f, g):
        """Polynomial pseudo-division of `f` and `g`. """
        dom, per, F, G = f.unify(g)
        q, r = dup_pdiv(F, G, dom)
        return per(q), per(r)

    def prem(f, g):
        """Polynomial pseudo-remainder of `f` and `g`. """
        dom, per, F, G = f.unify(g)
        return per(dup_prem(F, G, dom))

    def pquo(f, g):
        """Polynomial pseudo-quotient of `f` and `g`. """
        dom, per, F, G = f.unify(g)
        return per(dup_pquo(F, G, dom))

    def pexquo(f, g):
        """Polynomial exact pseudo-quotient of `f` and `g`. """
        dom, per, F, G = f.unify(g)
        return per(dup_pexquo(F, G, dom))

    def div(f, g):
        """Polynomial division with remainder of `f` and `g`. """
        dom, per, F, G = f.unify(g)
        q, r = dup_div(F, G, dom)
        return per(q), per(r)

    def rem(f, g):
        """Computes polynomial remainder in of `f` and `g`. """
        dom, per, F, G = f.unify(g)
        return per(dup_rem(F, G, dom))

    def quo(f, g):
        """Computes polynomial quotient in of `f` and `g`. """
        dom, per, F, G = f.unify(g)
        return per(dup_quo(F, G, dom))

    def exquo(f, g):
        """Computes polynomial exact quotient in of `f` and `g`. """
        dom, per, F, G = f.unify(g)
        return per(dup_exquo(F, G, dom))

    def degree(f):
        """Returns the leading degree of `f`. """
        return dup_degree(f.rep)

    def LC(f):
        """Returns the leading coefficent of `f`. """
        return dup_LC(f.rep, f.dom)

    def TC(f):
        """Returns the trailing coefficent of `f`. """
        return dup_TC(f.rep, f.dom)

    def nth(f, n):
        """Returns the `n`-th coefficient of `f`. """
        if isinstance(n, int):
            return dup_nth(f.rep, n, f.dom)
        else:
            raise TypeError("`int` expected, got %s" % type(n))

    def max_norm(f):
        """Returns maximum norm of `f`. """
        return dup_max_norm(f.rep, f.dom)

    def l1_norm(f):
        """Returns l1 norm of `f`. """
        return dup_l1_norm(f.rep, f.dom)

    def ground_to_ring(f):
        """Clear denominators, but keep the ground domain. """
        coeff, F = dup_ground_to_ring(f.rep, f.dom, f.dom.get_ring())
        return coeff, f.per(F)

    def integrate(f, m=1):
        """Computes indefinite integral of `f`. """
        if isinstance(m, int):
            return f.per(dup_integrate(f.rep, m, f.dom))
        else:
            raise TypeError("`int` expected, got %s" % type(m))

    def diff(f, m=1):
        """Computes partial derivative of `f`. """
        if isinstance(m, int):
            return f.per(dup_diff(f.rep, m, f.dom))
        else:
            raise TypeError("`int` expected, got %s" % type(m))

    def eval(f, a):
        """Evaluates `f` at the given point `a`. """
        return dup_eval(f.rep, f.dom.convert(a), f.dom)

    def half_gcdex(f, g):
        """Half extended Euclidean algorithm. """
        dom, per, F, G = f.unify(g)
        s, h = dup_half_gcdex(F, G, dom)
        return per(s), per(h)

    def gcdex(f, g):
        """Extended Euclidean algorithm. """
        dom, per, F, G = f.unify(g)
        s, t, h = dup_gcdex(F, G, dom)
        return per(s), per(t), per(h)

    def invert(f, g):
        """Invert `f` modulo `g`, if possible. """
        dom, per, F, G = f.unify(g)
        return per(dup_invert(F, G, dom))

    def subresultants(f, g, **flags):
        """Computes subresultant PRS sequence of `f` and `g`. """
        dom, per, F, G = f.unify(g)
        R = dup_subresultants(F, G, dom)
        return map(per, R)

    def resultant(f, g):
        """Computes resultant of `f` and `g` via PRS. """
        dom, _, F, G = f.unify(g)
        return dup_resultant(F, G, dom)

    def discriminant(f):
        """Computes discriminant of `f`. """
        return dup_discriminant(f.rep, f.dom)

    def cofactors(f, g):
        """Returns GCD of `f` and `g` and their cofactors. """
        dom, per, F, G = f.unify(g)
        h, cff, cfg = dup_inner_gcd(F, G, dom)
        return per(h), per(cff), per(cfg)

    def gcd(f, g):
        """Returns polynomial GCD of `f` and `g`. """
        dom, per, F, G = f.unify(g)
        return per(dup_gcd(F, G, dom))

    def lcm(f, g):
        """Returns polynomial LCM of `f` and `g`. """
        dom, per, F, G = f.unify(g)
        return per(dup_lcm(F, G, dom))

    def trunc(f, p):
        """Reduce `f` modulo a constant `p`. """
        return f.per(dup_trunc(f.rep, f.dom.convert(p), f.dom))

    def monic(f):
        """Divides all coefficients by `LC(f)`. """
        return f.per(dup_monic(f.rep, f.dom))

    def content(f):
        """Returns GCD of polynomial coefficients. """
        return dup_content(f.rep, f.dom)

    def primitive(f):
        """Returns content and a primitive form of `f`. """
        cont, F = dup_primitive(f.rep, f.dom)
        return cont, f.per(F)

    def compose(f, g):
        """Computes functional composition of `f` and `g`. """
        dom, per, F, G = f.unify(g)
        return per(dup_compose(F, G, dom))

    def decompose(f):
        """Computes functional decomposition of `f`. """
        return map(f.per, dup_decompose(f.rep, f.dom))

    def sturm(f):
        """Computes the Sturm sequence of `f`. """
        return map(f.per, dup_sturm(f.rep, f.dom))

    def sqf_part(f):
        """Computes square-free part of `f`. """
        return f.per(dup_sqf_part(f.rep, f.dom))

    def sqf_list(f, all=False):
        """Returns a list of square-free factors of `f`. """
        coeff, factors = dup_sqf_list(f.rep, f.dom, all)
        return coeff, [ (f.per(g), k) for g, k in factors ]

    def sqf_list_include(f, all=False):
        """Returns a list of square-free factors of `f`. """
        factors = dup_sqf_list_include(f.rep, f.dom, all)
        return [ (f.per(g), k) for g, k in factors ]

    def factor_list(f):
        """Returns a list of irreducible factors of `f`. """
        coeff, factors = dup_factor_list(f.rep, f.dom)
        return coeff, [ (f.per(g), k) for g, k in factors ]

    def factor_list_include(f):
        """Returns a list of irreducible factors of `f`. """
        factors = dup_factor_list_include(f.rep, f.dom)
        return [ (f.per(g), k) for g, k in factors ]

    @property
    def is_zero(f):
        """Returns `True` if `f` is a zero polynomial. """
        return not f

    @property
    def is_one(f):
        """Returns `True` if `f` is a unit polynomial. """
        return f.rep == [f.dom.one]

    @property
    def is_ground(f):
        """Returns `True` if `f` is an element of the ground domain. """
        return len(f.rep) <= 1

    @property
    def is_sqf(f):
        """Returns `True` if `f` is a square-free polynomial. """
        return dup_sqf_p(f.rep, f.dom)

    @property
    def is_monic(f):
        """Returns `True` if the leading coefficient of `f` is one. """
        return f.dom.is_one(dup_LC(f.rep, f.dom))

    @property
    def is_primitive(f):
        """Returns `True` if GCD of coefficients of `f` is one. """
        return f.dom.is_one(dup_content(f.rep, f.dom))

    def __abs__(f):
        return f.abs()

    def __neg__(f):
        return f.neg()

    def __add__(f, g):
        if isinstance(g, DUP):
            return f.add(g)
        else:
            try:
                return f.add_term(g, 0)
            except TypeError:
                return NotImplemented

    def __radd__(f, g):
        return f.__add__(g)

    def __sub__(f, g):
        if isinstance(g, DUP):
            return f.sub(g)
        else:
            try:
                return f.sub_term(g, 0)
            except TypeError:
                return NotImplemented

    def __rsub__(f, g):
        return (-f).__add__(g)

    def __mul__(f, g):
        if isinstance(g, DUP):
            return f.mul(g)
        else:
            try:
                return f.mul_ground(g)
            except TypeError:
                return NotImplemented

    def __rmul__(f, g):
        return f.__mul__(g)

    def __pow__(f, n):
        return f.pow(n)

    def __divmod__(f, g):
        return f.div(g)

    def __mod__(f, g):
        return f.rem(g)

    def __floordiv__(f, g):
        if isinstance(g, DUP):
            return f.exquo(g)
        else:
            try:
                return f.exquo_ground(g)
            except TypeError:
                return NotImplemented

    def __eq__(f, g):
        try:
            _, _, F, G = f.unify(g)

            return F == G
        except UnificationFailed:
            return False

    def __ne__(f, g):
        try:
            _, _, F, G = f.unify(g)

            return F != G
        except UnificationFailed:
            return True

    def __nonzero__(f):
        return bool(f.rep)

def init_normal_DMP(rep, lev, dom):
    return DMP(dmp_normal(rep, lev, dom), dom, lev)

class DMP(object):
    """Dense Multivariate Polynomials over `K`. """

    __slots__ = ['rep', 'lev', 'dom']

    def __init__(self, rep, dom, lev=None):
        if lev is not None:
            if type(rep) is dict:
                rep = dmp_from_dict(rep, lev, dom)
            elif type(rep) is not list:
                rep = dmp_ground(dom.convert(rep), lev)
        else:
            rep, lev = dmp_validate(rep)

        self.rep = rep
        self.lev = lev
        self.dom = dom

    def __repr__(f):
        return "%s(%s, %s)" % (f.__class__.__name__, f.rep, f.dom)

    def __hash__(f):
        return hash((f.__class__.__name__, repr(f.rep), f.lev, f.dom))

    def __getstate__(self):
        return (self.rep, self.lev, self.dom)

    def __getnewargs__(self):
        return (self.rep, self.lev, self.dom)

    def unify(f, g):
        """Unify representations of two multivariate polynomials. """
        if not isinstance(g, DMP) or f.lev != g.lev:
            raise UnificationFailed("can't unify %s with %s" % (f, g))

        if f.dom == g.dom:
            return f.lev, f.dom, f.per, f.rep, g.rep
        else:
            lev, dom = f.lev, f.dom.unify(g.dom)

            F = dmp_convert(f.rep, lev, f.dom, dom)
            G = dmp_convert(g.rep, lev, g.dom, dom)

            def per(rep, dom=dom, lev=lev, kill=False):
                if kill:
                    if not lev:
                        return rep
                    else:
                        lev -= 1

                return DMP(rep, dom, lev)

            return lev, dom, per, F, G

    def per(f, rep, dom=None, kill=False):
        """Create a DMP out of the given representation. """
        lev = f.lev

        if kill:
            if not lev:
                return rep
            else:
                lev -= 1

        if dom is None:
            dom = f.dom

        return DMP(rep, dom, lev)

    @classmethod
    def zero(cls, lev, dom):
        return DMP(0, dom, lev)

    @classmethod
    def one(cls, lev, dom):
        return DMP(1, dom, lev)

    @classmethod
    def from_list(cls, rep, lev, dom):
        """Create an instance of `cls` given a list of native coefficients. """
        return cls(dmp_convert(rep, lev, None, dom), dom, lev)

    @classmethod
    def from_sympy_list(cls, rep, lev, dom):
        """Create an instance of `cls` given a list of SymPy coefficients. """
        return cls(dmp_from_sympy(rep, lev, dom), dom, lev)

    def to_dict(f):
        """Convert `f` to a dict representation with native coefficients. """
        return dmp_to_dict(f.rep, f.lev)

    def to_sympy_dict(f):
        """Convert `f` to a dict representation with SymPy coefficients. """
        rep = dmp_to_dict(f.rep, f.lev)

        for k, v in rep.iteritems():
            rep[k] = f.dom.to_sympy(v)

        return rep

    def to_ring(f):
        """Make the ground domain a field. """
        return f.convert(f.dom.get_ring())

    def to_field(f):
        """Make the ground domain a field. """
        return f.convert(f.dom.get_field())

    def convert(f, dom):
        """Convert the ground domain of `f`. """
        if f.dom == dom:
            return f
        else:
            return DMP(dmp_convert(f.rep, f.lev, f.dom, dom), dom, f.lev)

    def coeffs(f):
        """Returns all non-zero coefficients from `f` in lex order. """
        return [ c for _, c in dmp_list_terms(f.rep, f.lev, f.dom) ]

    def monoms(f):
        """Returns all non-zero monomials from `f` in lex order. """
        return [ m for m, _ in dmp_list_terms(f.rep, f.lev, f.dom) ]

    def terms(f):
        """Returns all non-zero terms from `f` in lex order. """
        return dmp_list_terms(f.rep, f.lev, f.dom)

    def all_coeffs(f):
        """Returns all coefficients from `f`. """
        if not f.lev:
            if not f:
                return [f.dom.zero]
            else:
                return [ c for c in f.rep ]
        else:
            raise PolynomialError('multivariate polynomials not supported')

    def all_monoms(f):
        """Returns all monomials from `f`. """
        if not f.lev:
            n = gf_degree(f.rep)

            if n < 0:
                return [(0,)]
            else:
                return [ (n-i,) for i, c in enumerate(f.rep) ]
        else:
            raise PolynomialError('multivariate polynomials not supported')

    def all_terms(f):
        """Returns all terms from a `f`. """
        if not f.lev:
            n = gf_degree(f.rep)

            if n < 0:
                return [((0,), f.dom.zero)]
            else:
                return [ ((n-i,), c) for i, c in enumerate(f.rep) ]
        else:
            raise PolynomialError('multivariate polynomials not supported')

    def lift(f):
        """Convert algebraic coefficients to rationals. """
        return f.per(dmp_lift(f.rep, f.lev, f.dom), dom=f.dom.dom)

    def deflate(f):
        """Reduce degree of `f` by mapping `x_i**m` to `y_i`. """
        J, F = dmp_deflate(f.rep, f.lev, f.dom)
        return J, f.per(F)

    def terms_gcd(f):
        """Remove GCD of terms from the polynomial `f`. """
        J, F = dmp_terms_gcd(f.rep, f.lev, f.dom)
        return J, f.per(F)

    def mul_ground(f, c):
        """Multiply `f` by a an element of the ground domain. """
        return f.per(dmp_mul_ground(f.rep, f.dom.convert(c), f.lev, f.dom))

    def quo_ground(f, c):
        """Quotient of `f` by a an element of the ground domain. """
        return f.per(dmp_quo_ground(f.rep, f.dom.convert(c), f.lev, f.dom))

    def exquo_ground(f, c):
        """Exact quotient of `f` by a an element of the ground domain. """
        return f.per(dmp_exquo_ground(f.rep, f.dom.convert(c), f.lev, f.dom))

    def abs(f):
        """Make all coefficients in `f` positive. """
        return f.per(dmp_abs(f.rep, f.lev, f.dom))

    def neg(f):
        """Negate all cefficients in `f`. """
        return f.per(dmp_neg(f.rep, f.lev, f.dom))

    def add(f, g):
        """Add two multivariate polynomials `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_add(F, G, lev, dom))

    def sub(f, g):
        """Subtract two multivariate polynomials `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_sub(F, G, lev, dom))

    def mul(f, g):
        """Multiply two multivariate polynomials `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_mul(F, G, lev, dom))

    def sqr(f):
        """Square a multivariate polynomial `f`. """
        return f.per(dmp_sqr(f.rep, f.lev, f.dom))

    def pow(f, n):
        """Raise `f` to a non-negative power `n`. """
        if isinstance(n, int):
            return f.per(dmp_pow(f.rep, n, f.lev, f.dom))
        else:
            raise TypeError("`int` expected, got %s" % type(n))

    def pdiv(f, g):
        """Polynomial pseudo-division of `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        q, r = dmp_pdiv(F, G, lev, dom)
        return per(q), per(r)

    def prem(f, g):
        """Polynomial pseudo-remainder of `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_prem(F, G, lev, dom))

    def pquo(f, g):
        """Polynomial pseudo-quotient of `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_pquo(F, G, lev, dom))

    def pexquo(f, g):
        """Polynomial exact pseudo-quotient of `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_pexquo(F, G, lev, dom))

    def div(f, g):
        """Polynomial division with remainder of `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        q, r = dmp_div(F, G, lev, dom)
        return per(q), per(r)

    def rem(f, g):
        """Computes polynomial remainder of `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_rem(F, G, lev, dom))

    def quo(f, g):
        """Computes polynomial quotient of `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_quo(F, G, lev, dom))

    def exquo(f, g):
        """Computes polynomial exact quotient of `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_exquo(F, G, lev, dom))

    def degree(f, j=0):
        """Returns the leading degree of `f` in `x_j`. """
        if isinstance(j, int):
            return dmp_degree_in(f.rep, j, f.lev)
        else:
            raise TypeError("`int` expected, got %s" % type(j))

    def degree_list(f):
        """Returns a list of degrees of `f`. """
        return dmp_degree_list(f.rep, f.lev)

    def total_degree(f):
        """Returns the total degree of `f`. """
        return sum(dmp_degree_list(f.rep, f.lev))

    def LC(f):
        """Returns the leading coefficent of `f`. """
        return dmp_ground_LC(f.rep, f.lev, f.dom)

    def TC(f):
        """Returns the trailing coefficent of `f`. """
        return dmp_ground_TC(f.rep, f.lev, f.dom)

    def nth(f, *N):
        """Returns the `n`-th coefficient of `f`. """
        if all(isinstance(n, int) for n in N):
            return dmp_ground_nth(f.rep, N, f.lev, f.dom)
        else:
            raise TypeError("a sequence of integers expected")

    def max_norm(f):
        """Returns maximum norm of `f`. """
        return dmp_max_norm(f.rep, f.lev, f.dom)

    def l1_norm(f):
        """Returns l1 norm of `f`. """
        return dmp_l1_norm(f.rep, f.lev, f.dom)

    def ground_to_ring(f):
        """Clear denominators, but keep the ground domain. """
        coeff, F = dmp_ground_to_ring(f.rep, f.lev, f.dom, f.dom.get_ring())
        return coeff, f.per(F)

    def integrate(f, m=1, j=0):
        """Computes indefinite integral of `f`. """
        if not isinstance(m, int):
            raise TypeError("`int` expected, got %s" % type(m))

        if not isinstance(j, int):
            raise TypeError("`int` expected, got %s" % type(j))

        return f.per(dmp_integrate_in(f.rep, m, j, f.lev, f.dom))

    def diff(f, m=1, j=0):
        """Computes `m`-th order derivative of `f` in `x_j`. """
        if not isinstance(m, int):
            raise TypeError("`int` expected, got %s" % type(m))

        if not isinstance(j, int):
            raise TypeError("`int` expected, got %s" % type(j))

        return f.per(dmp_diff_in(f.rep, m, j, f.lev, f.dom))

    def eval(f, a, j=0):
        """Evaluates `f` at the given point `a` in `x_j`. """
        if not isinstance(j, int):
            raise TypeError("`int` expected, got %s" % type(j))

        return f.per(dmp_eval_in(f.rep,
            f.dom.convert(a), j, f.lev, f.dom), kill=True)

    def half_gcdex(f, g):
        """Half extended Euclidean algorithm, if univariate. """
        lev, dom, per, F, G = f.unify(g)

        if not lev:
            s, h = dup_half_gcdex(F, G, dom)
            return per(s), per(h)
        else:
            raise ValueError('univariate polynomial expected')

    def gcdex(f, g):
        """Extended Euclidean algorithm, if univariate. """
        lev, dom, per, F, G = f.unify(g)

        if not lev:
            s, t, h = dup_gcdex(F, G, dom)
            return per(s), per(t), per(h)
        else:
            raise ValueError('univariate polynomial expected')

    def invert(f, g):
        """Invert `f` modulo `g`, if possible. """
        lev, dom, per, F, G = f.unify(g)

        if not lev:
            return per(dup_invert(F, G, dom))
        else:
            raise ValueError('univariate polynomial expected')

    def subresultants(f, g):
        """Computes subresultant PRS sequence of `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        R = dmp_subresultants(F, G, lev, dom)
        return map(per, R)

    def resultant(f, g):
        """Computes resultant of `f` and `g` via PRS. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_resultant(F, G, lev, dom), kill=True)

    def discriminant(f):
        """Computes discriminant of `f`. """
        return f.per(dmp_discriminant(f.rep, f.lev, f.dom), kill=True)

    def cofactors(f, g):
        """Returns GCD of `f` and `g` and their cofactors. """
        lev, dom, per, F, G = f.unify(g)
        h, cff, cfg = dmp_inner_gcd(F, G, lev, dom)
        return per(h), per(cff), per(cfg)

    def gcd(f, g):
        """Returns polynomial GCD of `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_gcd(F, G, lev, dom))

    def lcm(f, g):
        """Returns polynomial LCM of `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_lcm(F, G, lev, dom))

    def trunc(f, p):
        """Reduce `f` modulo a constant `p`. """
        return f.per(dmp_ground_trunc(f.rep, f.dom.convert(p), f.lev, f.dom))

    def monic(f):
        """Divides all coefficients by `LC(f)`. """
        return f.per(dmp_ground_monic(f.rep, f.lev, f.dom))

    def content(f):
        """Returns GCD of polynomial coefficients. """
        return dmp_ground_content(f.rep, f.lev, f.dom)

    def primitive(f):
        """Returns content and a primitive form of `f`. """
        cont, F = dmp_ground_primitive(f.rep, f.lev, f.dom)
        return cont, f.per(F)

    def compose(f, g):
        """Computes functional composition of `f` and `g`. """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_compose(F, G, lev, dom))

    def decompose(f):
        """Computes functional decomposition of `f`. """
        if not f.lev:
            return map(f.per, dup_decompose(f.rep, f.dom))
        else:
            raise ValueError('univariate polynomial expected')

    def sturm(f):
        """Computes the Sturm sequence of `f`. """
        if not f.lev:
            return map(f.per, dup_sturm(f.rep, f.dom))
        else:
            raise ValueError('univariate polynomial expected')

    def sqf_norm(f):
        """Computes square-free norm of `f`. """
        s, g, r = dmp_sqf_norm(f.rep, f.lev, f.dom)
        return s, f.per(g), f.per(r, dom=f.dom.dom)

    def sqf_part(f):
        """Computes square-free part of `f`. """
        return f.per(dmp_sqf_part(f.rep, f.lev, f.dom))

    def sqf_list(f, all=False):
        """Returns a list of square-free factors of `f`. """
        coeff, factors = dmp_sqf_list(f.rep, f.lev, f.dom, all)
        return coeff, [ (f.per(g), k) for g, k in factors ]

    def sqf_list_include(f, all=False):
        """Returns a list of square-free factors of `f`. """
        factors = dmp_sqf_list_include(f.rep, f.lev, f.dom, all)
        return [ (f.per(g), k) for g, k in factors ]

    def factor_list(f):
        """Returns a list of irreducible factors of `f`. """
        coeff, factors = dmp_factor_list(f.rep, f.lev, f.dom)
        return coeff, [ (f.per(g), k) for g, k in factors ]

    def factor_list_include(f):
        """Returns a list of irreducible factors of `f`. """
        factors = dmp_factor_list_include(f.rep, f.lev, f.dom)
        return [ (f.per(g), k) for g, k in factors ]

    def intervals_sqf(f, eps=None, fast=False):
        """Compute isolating intervals for roots of square-free `f`. """
        if not f.lev:
            return dup_isolate_real_roots_sqf(f.rep, f.dom, eps=eps, fast=fast)
        else:
            raise PolynomialError("can't isolate roots of a multivariate polynomial")

    def intervals(f, eps=None, fast=False):
        """Compute isolating intervals for roots of `f`. """
        if not f.lev:
            return dup_isolate_real_roots(f.rep, f.dom, eps=eps, fast=fast)
        else:
            raise PolynomialError("can't isolate roots of a multivariate polynomial")

    @property
    def is_zero(f):
        """Returns `True` if `f` is a zero polynomial. """
        return dmp_zero_p(f.rep, f.lev)

    @property
    def is_one(f):
        """Returns `True` if `f` is a unit polynomial. """
        return dmp_one_p(f.rep, f.lev, f.dom)

    @property
    def is_ground(f):
        """Returns `True` if `f` is an element of the ground domain. """
        return dmp_ground_p(f.rep, None, f.lev)

    @property
    def is_sqf(f):
        """Returns `True` if `f` is a square-free polynomial. """
        return dmp_sqf_p(f.rep, f.lev, f.dom)

    @property
    def is_monic(f):
        """Returns `True` if the leading coefficient of `f` is one. """
        return f.dom.is_one(dmp_ground_LC(f.rep, f.lev, f.dom))

    @property
    def is_primitive(f):
        """Returns `True` if GCD of coefficients of `f` is one. """
        return f.dom.is_one(dmp_ground_content(f.rep, f.lev, f.dom))

    @property
    def is_ground(f):
        """Returns `True` if `f` is an element of the ground domain. """
        return all(d <= 0 for d in dmp_degree_list(f.rep, f.lev))

    @property
    def is_linear(f):
        """Returns `True` if `f` is linear in all its variables. """
        return all(d <= 1 for d in dmp_degree_list(f.rep, f.lev))

    @property
    def is_homogeneous(f):
        """Returns `True` if `f` has zero trailing coefficient. """
        return f.dom.is_zero(dmp_ground_TC(f.rep, f.lev, f.dom))

    def __abs__(f):
        return f.abs()

    def __neg__(f):
        return f.neg()

    def __add__(f, g):
        if not isinstance(g, DMP):
            try:
                g = f.per(dmp_ground(f.dom.convert(g), f.lev))
            except TypeError:
                return NotImplemented

        return f.add(g)

    def __radd__(f, g):
        return f.__add__(g)

    def __sub__(f, g):
        if not isinstance(g, DMP):
            try:
                g = f.per(dmp_ground(f.dom.convert(g), f.lev))
            except TypeError:
                return NotImplemented

        return f.sub(g)

    def __rsub__(f, g):
        return (-f).__add__(g)

    def __mul__(f, g):
        if isinstance(g, DMP):
            return f.mul(g)
        else:
            try:
                return f.mul_ground(g)
            except TypeError:
                return NotImplemented

    def __rmul__(f, g):
        return f.__mul__(g)

    def __pow__(f, n):
        return f.pow(n)

    def __divmod__(f, g):
        return f.div(g)

    def __mod__(f, g):
        return f.rem(g)

    def __floordiv__(f, g):
        if isinstance(g, DMP):
            return f.exquo(g)
        else:
            try:
                return f.exquo_ground(g)
            except TypeError:
                return NotImplemented

    def __eq__(f, g):
        try:
            _, _, _, F, G = f.unify(g)

            if f.lev == g.lev:
                return F == G
        except UnificationFailed:
            pass

        return False

    def __ne__(f, g):
        try:
            _, _, _, F, G = f.unify(g)

            if f.lev == g.lev:
                return F != G
        except UnificationFailed:
            pass

        return True

    def __nonzero__(f):
        return not dmp_zero_p(f.rep, f.lev)

def init_normal_SDP(rep, lev, order, dom):
    raise NotImplementedError

class SDP(object):
    """Sparse Distributed Polynomials over `K`. """

    __slots__ = ['rep', 'lev', 'order', 'dom']

    # TODO: first make sdp_* functions more useful

def init_normal_DMF(num, den, lev, dom):
    return DFP(dmp_normal(num, lev, dom),
               dmp_normal(den, lev, dom), dom, lev)

class DMF(object):
    """Dense Multivariate Fractions over `K`. """

    __slots__ = ['num', 'den', 'lev', 'dom']

    def __init__(self, rep, dom, lev=None):
        assert dom.has_Ring, "QQ in ground not supported, yet"

        if type(rep) is tuple:
            num, den = rep

            if lev is not None:
                if type(num) is dict:
                    num = dmp_from_dict(num, lev, dom)

                if type(den) is dict:
                    den = dmp_from_dict(den, lev, dom)
            else:
                num, num_lev = dmp_validate(num)
                den, den_lev = dmp_validate(den)

                if num_lev == den_lev:
                    lev = num_lev
                else:
                    raise ValueError('inconsistent number of levels')

            if dmp_zero_p(den, lev):
                raise ZeroDivisionError('fraction denominator')

            if dmp_zero_p(num, lev):
                den = dmp_one(lev, dom)
            else:
                if dmp_negative_p(den, lev, dom):
                    num = dmp_neg(num, lev, dom)
                    den = dmp_neg(den, lev, dom)
        else:
            num = rep

            if lev is not None:
                if type(num) is dict:
                    num = dmp_from_dict(num, lev, dom)
                elif type(num) is not list:
                    num = dmp_ground(dom.convert(num), lev)
            else:
                num, lev = dmp_validate(num)

            den = dmp_one(lev, dom)

        self.num = num
        self.den = den
        self.lev = lev
        self.dom = dom

    def __repr__(f):
        return "%s((%s, %s), %s)" % (f.__class__.__name__, f.num, f.den, f.dom)

    def __hash__(f):
        return hash((f.__class__.__name__, repr(f.num), repr(f.den), f.lev, f.dom))

    def __getstate__(self):
        return (self.num, self.den, self.lev, self.dom)

    def __getnewargs__(self):
        return (self.num, self.den, self.lev, self.dom)

    def poly_unify(f, g):
        """Unify a multivariate fraction and a polynomial. """
        if not isinstance(g, DMP) or f.lev != g.lev:
            raise UnificationFailed("can't unify %s with %s" % (f, g))

        if f.dom == g.dom:
            return (f.lev, f.dom, f.per, (f.num, f.den), g.rep)
        else:
            lev, dom = f.lev, f.dom.unify(g.dom)

            F = (dmp_convert(f.num, lev, f.dom, dom),
                 dmp_convert(f.den, lev, f.dom, dom))

            G = dmp_convert(g.rep, lev, g.dom, dom)

            def per(num, den, dom=dom, lev=lev,
                    cancel=True, kill=False):
                if kill:
                    if not lev:
                        return num/den
                    else:
                        lev -= 1

                if cancel:
                    _, num, den = dmp_inner_gcd(num, den, lev, dom)

                return DMF((num, den), dom, lev)

            return lev, dom, per, F, G

    def frac_unify(f, g):
        """Unify representations of two multivariate fractions. """
        if not isinstance(g, DMF) or f.lev != g.lev:
            raise UnificationFailed("can't unify %s with %s" % (f, g))

        if f.dom == g.dom:
            return (f.lev, f.dom, f.per, (f.num, f.den),
                                         (g.num, g.den))
        else:
            lev, dom = f.lev, f.dom.unify(g.dom)

            F = (dmp_convert(f.num, lev, f.dom, dom),
                 dmp_convert(f.den, lev, f.dom, dom))

            G = (dmp_convert(g.num, lev, g.dom, dom),
                 dmp_convert(g.den, lev, g.dom, dom))

            def per(num, den, dom=dom, lev=lev,
                    cancel=True, kill=False):
                if kill:
                    if not lev:
                        return num/den
                    else:
                        lev -= 1

                if cancel:
                    _, num, den = dmp_inner_gcd(num, den, lev, dom)

                return DMF((num, den), dom, lev)

            return lev, dom, per, F, G

    def per(f, num, den, cancel=True, kill=False):
        """Create a DMF out of the given representation. """
        lev, dom = f.lev, f.dom

        if kill:
            if not lev:
                return num/den
            else:
                lev -= 1

        if cancel:
            _, num, den = dmp_inner_gcd(num, den, lev, dom)

        return DMF((num, den), dom, lev)

    def half_per(f, rep, kill=False):
        """Create a DMP out of the given representation. """
        lev = f.lev

        if kill:
            if not lev:
                return rep
            else:
                lev -= 1

        return DMP(rep, f.dom, lev)

    @classmethod
    def zero(cls, lev, dom):
        return DMF(0, dom, lev)

    @classmethod
    def one(cls, lev, dom):
        return DMF(1, dom, lev)

    def numer(f):
        """Returns numerator of `f`. """
        return f.half_per(f.num)

    def denom(f):
        """Returns denominator of `f`. """
        return f.half_per(f.den)

    def cancel(f):
        """Remove common factors from `f.num` and `f.den`. """
        return f.per(f.num, f.den)

    def neg(f):
        """Negate all cefficients in `f`. """
        return f.per(dmp_neg(f.num, f.lev, f.dom), f.den, cancel=False)

    def add(f, g):
        """Add two multivariate fractions `f` and `g`. """
        if isinstance(g, DMP):
            lev, dom, per, (F_num, F_den), G = f.poly_unify(g)
            num, den = dmp_add_mul(F_num, F_den, G, lev, dom), F_den
        else:
            lev, dom, per, F, G = f.frac_unify(g)
            (F_num, F_den), (G_num, G_den) = F, G

            num = dmp_add(dmp_mul(F_num, G_den, lev, dom),
                          dmp_mul(F_den, G_num, lev, dom), lev, dom)
            den = dmp_mul(F_den, G_den, lev, dom)

        return per(num, den)

    def sub(f, g):
        """Subtract two multivariate fractions `f` and `g`. """
        if isinstance(g, DMP):
            lev, dom, per, (F_num, F_den), G = f.poly_unify(g)
            num, den = dmp_sub_mul(F_num, F_den, G, lev, dom), F_den
        else:
            lev, dom, per, F, G = f.frac_unify(g)
            (F_num, F_den), (G_num, G_den) = F, G

            num = dmp_sub(dmp_mul(F_num, G_den, lev, dom),
                          dmp_mul(F_den, G_num, lev, dom), lev, dom)
            den = dmp_mul(F_den, G_den, lev, dom)

        return per(num, den)

    def mul(f, g):
        """Multiply two multivariate fractions `f` and `g`. """
        if isinstance(g, DMP):
            lev, dom, per, (F_num, F_den), G = f.poly_unify(g)
            num, den = dmp_mul(F_num, G, lev, dom), F_den
        else:
            lev, dom, per, F, G = f.frac_unify(g)
            (F_num, F_den), (G_num, G_den) = F, G

            num = dmp_mul(F_num, G_num, lev, dom)
            den = dmp_mul(F_den, G_den, lev, dom)

        return per(num, den)

    def pow(f, n):
        """Raise `f` to a non-negative power `n`. """
        if isinstance(n, int):
            return f.per(dmp_pow(f.num, n, f.lev, f.dom),
                         dmp_pow(f.den, n, f.lev, f.dom), cancel=False)
        else:
            raise TypeError("`int` expected, got %s" % type(n))

    def quo(f, g):
        """Computes quotient of fractions `f` and `g`. """
        if isinstance(g, DMP):
            lev, dom, per, (F_num, F_den), G = f.poly_unify(g)
            num, den = F_num, dmp_mul(F_den, G, lev, dom)
        else:
            lev, dom, per, F, G = f.frac_unify(g)
            (F_num, F_den), (G_num, G_den) = F, G

            num = dmp_mul(F_num, G_den, lev, dom)
            den = dmp_mul(F_den, G_num, lev, dom)

        return per(num, den)

    exquo = quo

    def invert(f):
        """Computes inverse of a fraction `f`. """
        return f.per(f.den, f.num, cancel=False)

    @property
    def is_zero(f):
        """Returns `True` if `f` is a zero fraction. """
        return dmp_zero_p(f.num, f.lev)

    @property
    def is_one(f):
        """Returns `True` if `f` is a unit fraction. """
        return dmp_one_p(f.num, f.lev, f.dom) and \
               dmp_one_p(f.den, f.lev, f.dom)

    def __neg__(f):
        return f.neg()

    def __add__(f, g):
        if isinstance(g, (DMP, DMF)):
            return f.add(g)

        try:
            return f.add(f.half_per(g))
        except TypeError:
            return NotImplemented

    def __radd__(f, g):
        return f.__add__(g)

    def __sub__(f, g):
        if isinstance(g, (DMP, DMF)):
            return f.sub(g)

        try:
            return f.sub(f.half_per(g))
        except TypeError:
            return NotImplemented

    def __rsub__(f, g):
        return (-f).__add__(g)

    def __mul__(f, g):
        if isinstance(g, (DMP, DMF)):
            return f.mul(g)

        try:
            return f.mul(f.half_per(g))
        except TypeError:
            return NotImplemented

    def __rmul__(f, g):
        return f.__mul__(g)

    def __pow__(f, n):
        return f.pow(n)

    def __div__(f, g):
        if isinstance(g, (DMP, DMF)):
            return f.exquo(g)

        try:
            return f.exquo(f.half_per(g))
        except TypeError:
            return NotImplemented

    __truediv__ = __div__

    def __eq__(f, g):
        try:
            if isinstance(g, DMP):
                _, _, _, (F_num, F_den), G = f.poly_unify(g)

                if f.lev == g.lev:
                    return dmp_one_p(F_den, f.lev, f.dom) and F_num == G
            else:
                _, _, _, F, G = f.frac_unify(g)

                if f.lev == g.lev:
                    return F == G
        except UnificationFailed:
            pass

        return False

    def __ne__(f, g):
        try:
            if isinstance(g, DMP):
                _, _, _, (F_num, F_den), G = f.poly_unify(g)

                if f.lev == g.lev:
                    return not (dmp_one_p(F_den, f.lev, f.dom) and F_num == G)
            else:
                _, _, _, F, G = f.frac_unify(g)

                if f.lev == g.lev:
                    return F != G
        except UnificationFailed:
            pass

        return True

    def __nonzero__(f):
        return not dmp_zero_p(f.num, f.lev)

def init_normal_ANP(rep, mod, dom):
    return ANP(dup_normal(rep, dom),
               dup_normal(mod, dom), dom)

class ANP(object):
    """Dense Algebraic Number Polynomials over a field. """

    __slots__ = ['rep', 'mod', 'dom']

    def __init__(self, rep, mod, dom):
        if type(rep) is dict:
            self.rep = dup_from_dict(rep, dom)
        else:
            if type(rep) is not list:
                rep = [dom.convert(rep)]

            self.rep = dup_strip(rep)

        if isinstance(mod, (DUP, DMP)):
            self.mod = mod.rep
        else:
            if type(mod) is dict:
                self.mod = dup_from_dict(mod, dom)
            else:
                self.mod = dup_strip(mod)

        self.dom = dom

    def __repr__(f):
        return "%s(%s, %s, %s)" % (f.__class__.__name__, f.rep, f.mod, f.dom)

    def __hash__(f):
        return hash((f.__class__.__name__, repr(f.rep), f.mod, f.dom))

    def __getstate__(self):
        return (self.rep, self.mod, self.dom)

    def __getnewargs__(self):
        return (self.rep, self.mod, self.dom)

    def __cmp__(f, g):
        """Make sorting deterministic. """
        k = len(f.rep) - len(g.rep)

        if not k:
            return cmp(f.rep, g.rep)
        else:
            return k

    def unify(f, g):
        """Unify representations of two algebraic numbers. """
        if not isinstance(g, ANP) or f.mod != g.mod:
            raise UnificationFailed("can't unify %s with %s" % (f, g))

        if f.dom == g.dom:
            return f.dom, f.per, f.rep, g.rep, f.mod
        else:
            dom = f.dom.unify(g.dom)

            F = dup_convert(f.rep, f.dom, dom)
            G = dup_convert(g.rep, g.dom, dom)

            if dom != f.dom and dom != g.dom:
                mod = dup_convert(f.mod, f.dom, dom)
            else:
                if dom == f.dom:
                    H = f.mod
                else:
                    H = g.mod

            per = lambda rep: ANP(rep, mod, dom)

        return dom, per, F, G, mod

    def per(f, rep, mod=None, dom=None):
        return ANP(rep, mod or f.mod, dom or f.dom)

    @classmethod
    def zero(cls, mod, dom):
        return ANP(0, mod, dom)

    @classmethod
    def one(cls, mod, dom):
        return ANP(1, mod, dom)

    def to_dict(f):
        """Convert `f` to a dict representation with native coefficients. """
        return dmp_to_dict(f.rep, 0)

    def to_sympy_dict(f):
        """Convert `f` to a dict representation with SymPy coefficients. """
        rep = dmp_to_dict(f.rep, 0)

        for k, v in rep.iteritems():
            rep[k] = f.dom.to_sympy(v)

        return rep

    def to_list(f):
        """Convert `f` to a list representation with native coefficients. """
        return f.rep

    def to_sympy_list(f):
        """Convert `f` to a list representation with SymPy coefficients. """
        return [ f.dom.to_sympy(c) for c in f.rep ]

    @classmethod
    def from_list(cls, rep, mod, dom):
        return ANP(dup_strip(map(dom.convert, rep)), mod, dom)

    def neg(f):
        return f.per(dup_neg(f.rep, f.dom))

    def add(f, g):
        dom, per, F, G, mod = f.unify(g)
        return per(dup_add(F, G, dom))

    def sub(f, g):
        dom, per, F, G, mod = f.unify(g)
        return per(dup_sub(F, G, dom))

    def mul(f, g):
        dom, per, F, G, mod = f.unify(g)
        return per(dup_rem(dup_mul(F, G, dom), mod, dom))

    def pow(f, n):
        """Raise `f` to a non-negative power `n`. """
        if isinstance(n, int):
            if n < 0:
                F, n = dup_invert(f.rep, f.mod, f.dom), -n
            else:
                F = f.rep

            return f.per(dup_rem(dup_pow(F, n, f.dom), f.mod, f.dom))
        else:
            raise TypeError("`int` expected, got %s" % type(n))

    def div(f, g):
        dom, per, F, G, mod = f.unify(g)
        return (per(dup_rem(dup_mul(F, dup_invert(G, mod, dom), dom), mod, dom)), self.zero(mod, dom))

    def rem(f, g):
        dom, _, _, _, mod = f.unify(g)
        return self.zero(mod, dom)

    def quo(f, g):
        dom, per, F, G, mod = f.unify(g)
        return per(dup_rem(dup_mul(F, dup_invert(G, mod, dom), dom), mod, dom))

    exquo = quo

    def LC(f):
        """Returns the leading coefficent of `f`. """
        return dup_LC(f.rep, f.dom)

    def TC(f):
        """Returns the trailing coefficent of `f`. """
        return dup_TC(f.rep, f.dom)

    @property
    def is_zero(f):
        """Returns `True` if `f` is a zero algebraic number. """
        return not f

    @property
    def is_one(f):
        """Returns `True` if `f` is a unit algebraic number. """
        return f.rep == [f.dom.one]

    @property
    def is_ground(f):
        """Returns `True` if `f` is an element of the ground domain. """
        return not f.rep or len(f.rep) == 1

    def __neg__(f):
        return f.neg()

    def __add__(f, g):
        if isinstance(g, ANP):
            return f.add(g)
        else:
            try:
                return f.add(f.per(g))
            except TypeError:
                return NotImplemented

    def __radd__(f, g):
        return f.__add__(g)

    def __sub__(f, g):
        if isinstance(g, ANP):
            return f.sub(g)
        else:
            try:
                return f.sub(f.per(g))
            except TypeError:
                return NotImplemented

    def __rsub__(f, g):
        return (-f).__add__(g)

    def __mul__(f, g):
        if isinstance(g, ANP):
            return f.mul(g)
        else:
            try:
                return f.mul(f.per(g))
            except TypeError:
                return NotImplemented

    def __rmul__(f, g):
        return f.__mul__(g)

    def __pow__(f, n):
        return f.pow(n)

    def __divmod__(f, g):
        return f.div(g)

    def __mod__(f, g):
        return f.rem(g)

    def __div__(f, g):
        if isinstance(g, ANP):
            return f.exquo(g)
        else:
            try:
                return f.exquo(f.per(g))
            except TypeError:
                return NotImplemented

    __truediv__ = __div__

    def __eq__(f, g):
        try:
            _, _, F, G, _ = f.unify(g)

            return F == G
        except UnificationFailed:
            return False

    def __ne__(f, g):
        try:
            _, _, F, G, _ = f.unify(g)

            return F != G
        except UnificationFailed:
            return True

    def __nonzero__(f):
        return bool(f.rep)

