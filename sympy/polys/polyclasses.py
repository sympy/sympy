"""OO layer for several polynomial representations. """

class GenericPoly(object):
    """Base class for low--level polynomial representations. """

    def ground_to_ring(f):
        """Make the ground domain a ring. """
        return f.set_domain(f.dom.get_ring())

    def ground_to_field(f):
        """Make the ground domain a field. """
        return f.set_domain(f.dom.get_field())

    def ground_to_exact(f):
        """Make the ground domain exact. """
        return f.set_domain(f.dom.get_exact())

    @classmethod
    def _perify_factors(per, result, include):
        if include:
            coeff, factors = result
        else:
            coeff = result

        factors = [ (per(g), k) for g, k in factors ]

        if include:
            return coeff, factors
        else:
            return factors

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
    dup_clear_denoms, dmp_clear_denoms,
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
)

from sympy.polys.rootisolation import (
    dup_isolate_real_roots_sqf,
    dup_isolate_real_roots,
    dup_isolate_all_roots_sqf,
    dup_isolate_all_roots,
    dup_refine_real_root,
    dup_count_real_roots,
    dup_count_complex_roots,
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
    gf_gcdex, gf_gcd, gf_lcm, gf_cofactors,
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
        """
        Returns a zero polynomial with modulus `mod` and domain `dom`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP.zero(ZZ, 2)
        GFP([], 2, ZZ)
        """
        return GFP(0, mod, dom)

    @classmethod
    def one(cls, dom, mod):
        r"""
        Returns a one polynomial with modulus `mod` and domain `dom`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP.one(ZZ, 2) == \
        ... GFP([ZZ(1)], 2, ZZ)
        True
        """
        return GFP(1, mod, dom)

    def trunc(f, mod):
        r"""
        Reduce `f` using new modulus.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(7), ZZ(-2), ZZ(3)], 11, ZZ).trunc(5) == \
        ... GFP([ZZ(2), ZZ(4), ZZ(3)], 5, ZZ)
        True
        """
        if mod == f.mod:
            return f
        else:
            return GFP(gf_trunc(f.rep, mod), mod, f.dom, f.sym)

    def convert(f, dom):
        r"""
        Convert the ground domain of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(4), ZZ(0)], 3, ZZ).convert(ZZ) == \
        ... GFP([ZZ(1), ZZ(4), ZZ(0)], 3, ZZ)
        True
        """
        if f.dom == dom:
            return f
        elif dom.is_ZZ:
            return GFP(gf_convert(f.rep, f.mod, f.dom, dom), f.mod, dom, f.sym)
        else:
            raise DomainError("can't convert GFP ground domain to %s" % dom)

    def coeffs(f, order=None):
        """
        Returns all non-zero coefficients from `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(2), ZZ(0), ZZ(4)], 5, ZZ).coeffs()
        [1, 2, -1]
        """
        if not f:
            return [f.dom.zero]
        elif not f.sym:
            return [ c for c in f.rep if c ]
        else:
            return [ gf_int(c, f.mod) for c in f.rep if c ]

    def monoms(f, order=None):
        """
        Returns all non-zero monomials from `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(2), ZZ(0), ZZ(4)], 5, ZZ).monoms()
        [(3,), (2,), (0,)]
        """
        n = gf_degree(f.rep)

        if n < 0:
            return [(0,)]
        else:
            return [ (n-i,) for i, c in enumerate(f.rep) if c ]

    def terms(f, order=None):
        """
        Returns all non-zero terms from `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(2), ZZ(0), ZZ(4)], 5, ZZ).terms()
        [((3,), 1), ((2,), 2), ((0,), -1)]
        """
        n = gf_degree(f.rep)

        if n < 0:
            return [((0,), f.dom.zero)]
        elif not f.sym:
            return [ ((n-i,), c) for i, c in enumerate(f.rep) if c ]
        else:
            return [ ((n-i,), gf_int(c, f.mod)) for i, c in enumerate(f.rep) if c ]

    def all_coeffs(f):
        """
        Returns all coefficients from `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(2), ZZ(0), ZZ(4)], 5, ZZ).all_coeffs()
        [1, 2, 0, -1]
        """
        if not f:
            return [f.dom.zero]
        elif not f.sym:
            return [ c for c in f.rep ]
        else:
            return [ gf_int(c, f.mod) for c in f.rep ]

    def all_monoms(f):
        """
        Returns all monomials from `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(2), ZZ(0), ZZ(4)], 5, ZZ).all_monoms()
        [(3,), (2,), (1,), (0,)]
        """
        n = gf_degree(f.rep)

        if n < 0:
            return [((0,), f.dom.zero)]
        else:
            return [ (n-i,) for i, c in enumerate(f.rep) ]

    def all_terms(f):
        """
        Returns all terms from `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(2), ZZ(0), ZZ(4)], 5, ZZ).all_terms()
        [((3,), 1), ((2,), 2), ((1,), 0), ((0,), -1)]
        """
        n = gf_degree(f.rep)

        if n < 0:
            return [((0,), f.dom.zero)]
        elif not f.sym:
            return [ ((n-i,), c) for i, c in enumerate(f.rep) ]
        else:
            return [ ((n-i,), gf_int(c, f.mod)) for i, c in enumerate(f.rep) ]

    def deflate(f):
        r"""
        Reduce degree of `f` by mapping `x**m` to `y`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(0), ZZ(0), ZZ(1)], 5, ZZ).deflate() == \
        ... (3, GFP([ZZ(1), ZZ(1)], 5, ZZ))
        True
        """
        j, F = dup_deflate(f.rep, f.dom)
        return j, f.per(F)

    def terms_gcd(f):
        r"""
        Remove GCD of terms from the polynomial `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(0), ZZ(1), ZZ(0), ZZ(0)], 5, ZZ).terms_gcd() == \
        ... (2, GFP([ZZ(1), ZZ(0), ZZ(1)], 5, ZZ))
        True
        """
        j, F = dup_terms_gcd(f.rep, f.dom)
        return j, f.per(F)

    def add_ground(f, c):
        r"""
        Add an element of the ground domain to `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).add_ground(2) == \
        ... GFP([ZZ(3), ZZ(2), ZZ(1)], 5, ZZ)
        True
        """
        return f.per(gf_add_ground(f.rep, f.dom.convert(c), f.mod, f.dom))

    def sub_ground(f, c):
        r"""
        Subtract an element of the ground domain from `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).sub_ground(2) == \
        ... GFP([ZZ(3), ZZ(2), ZZ(2)], 5, ZZ)
        True
        """
        return f.per(gf_sub_ground(f.rep, f.dom.convert(c), f.mod, f.dom))

    def mul_ground(f, c):
        r"""
        Multiply `f` by an element of the ground domain.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).mul_ground(2) == \
        ... GFP([ZZ(1), ZZ(4), ZZ(3)], 5, ZZ)
        True
        """
        return f.per(gf_mul_ground(f.rep, f.dom.convert(c), f.mod, f.dom))

    def exquo_ground(f, c):
        r"""
        Divide `f` by an element of the ground domain.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).exquo_ground(2) == \
        ... GFP([ZZ(4), ZZ(1), ZZ(2)], 5, ZZ)
        True
        """
        return f.per(gf_exquo_ground(f.rep, f.dom.convert(c), f.mod, f.dom))

    def neg(f):
        r"""
        Negate all cefficients in `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(1), ZZ(0)], 5, ZZ).neg() == \
        ... GFP([ZZ(2), ZZ(3), ZZ(4), ZZ(0)], 5, ZZ)
        True
        """
        return f.per(gf_neg(f.rep, f.mod, f.dom))

    def add(f, g):
        r"""
        Add two univariate polynomials `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).add(
        ... GFP([ZZ(2), ZZ(2), ZZ(2)], 5, ZZ)) == \
        ... GFP([ZZ(4), ZZ(1)], 5, ZZ)
        True
        """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_add(F, G, mod, dom))

    def sub(f, g):
        r"""
        Subtract two univariate polynomials `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).sub(
        ... GFP([ZZ(2), ZZ(2), ZZ(2)], 5, ZZ)) == \
        ... GFP([ZZ(1), ZZ(0), ZZ(2)], 5, ZZ)
        True
        """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_sub(F, G, mod, dom))

    def mul(f, g):
        r"""
        Multiply two univariate polynomials `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).mul(
        ... GFP([ZZ(2), ZZ(2), ZZ(2)], 5, ZZ)) == \
        ... GFP([ZZ(1), ZZ(0), ZZ(3), ZZ(2), ZZ(3)], 5, ZZ)
        True
        """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_mul(F, G, mod, dom))

    def sqr(f):
        r"""
        Square a univariate polynomial `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).sqr() == \
        ... GFP([ZZ(4), ZZ(2), ZZ(3), ZZ(1), ZZ(1)], 5, ZZ)
        True
        """
        return f.per(gf_sqr(f.rep, f.mod, f.dom))

    def pow(f, n):
        r"""
        Raise `f` to a non-negative power `n`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).pow(3) == \
        ... GFP([ZZ(2), ZZ(4), ZZ(4), ZZ(2), ZZ(2), ZZ(1), ZZ(4)], 5, ZZ)
        True
        """
        if isinstance(n, int):
            return f.per(gf_pow(f.rep, n, f.mod, f.dom))
        else:
            raise TypeError("`int` expected, got %s" % type(n))

    def div(f, g):
        r"""
        Polynomial division with remainder of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(0), ZZ(1), ZZ(1)], 2, ZZ).div(
        ... GFP([ZZ(1), ZZ(1), ZZ(0)], 2, ZZ)) == \
        ... (GFP([ZZ(1), ZZ(1)], 2, ZZ), GFP([ZZ(1)], 2, ZZ))
        True
        """
        mod, dom, per, F, G = f.unify(g)
        q, r = gf_div(F, G, mod, dom)
        return per(q), per(r)

    def rem(f, g):
        r"""
        Computes polynomial remainder in of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(0), ZZ(1), ZZ(1)], 2, ZZ).rem(
        ... GFP([ZZ(1), ZZ(1), ZZ(0)], 2, ZZ)) == \
        ... GFP([ZZ(1)], 2, ZZ)
        True
        """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_rem(F, G, mod, dom))

    def quo(f, g):
        r"""
        Computes polynomial quotient in of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ

        >> GFP([ZZ(1), ZZ(0), ZZ(1), ZZ(1)], 2, ZZ).quo(
        ... GFP([ZZ(1), ZZ(1), ZZ(0)], 2, ZZ))
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: [1, 1, 0] does not divide [1, 0, 1, 1]
        >>> GFP([ZZ(1), ZZ(0), ZZ(3), ZZ(2), ZZ(3)], 5, ZZ).quo(
        ... GFP([ZZ(2), ZZ(2), ZZ(2)], 5, ZZ)) == \
        ... GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ)
        True
        """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_quo(F, G, mod, dom))

    def exquo(f, g):
        r"""
        Computes polynomial exact quotient in of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(0), ZZ(1), ZZ(1)], 2, ZZ).exquo(
        ... GFP([ZZ(1), ZZ(1), ZZ(0)], 2, ZZ)) == \
        ... GFP([ZZ(1), ZZ(1)], 2, ZZ)
        True
        >>> GFP([ZZ(1), ZZ(0), ZZ(3), ZZ(2), ZZ(3)], 5, ZZ).exquo(
        ... GFP([ZZ(2), ZZ(2), ZZ(2)], 5, ZZ)) == \
        ... GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ)
        True
        """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_exquo(F, G, mod, dom))

    def degree(f):
        """
        Returns the leading degree of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(1), ZZ(2), ZZ(0)], 5, ZZ).degree()
        3
        >>> GFP([], 5, ZZ).degree()
        -1
        """
        return gf_degree(f.rep)

    def LC(f):
        """
        Returns the leading coefficent of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(0), ZZ(1)], 5, ZZ).LC()
        3
        """
        return gf_LC(f.rep, f.dom)

    def TC(f):
        """
        Returns the trailing coefficent of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(0), ZZ(1)], 5, ZZ).TC()
        1
        """
        return gf_TC(f.rep, f.dom)

    def gcdex(f, g):
        r"""
        Extended Euclidean algorithm.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(8), ZZ(7)], 11, ZZ).gcdex(
        ... GFP([ZZ(1), ZZ(7), ZZ(1), ZZ(7)], 11, ZZ)) == \
        ... (GFP([ZZ(5), ZZ(6)], 11, ZZ), GFP([ZZ(6)], 11, ZZ),
        ...  GFP([ZZ(1), ZZ(7)], 11, ZZ))
        True
        """
        mod, dom, per, F, G = f.unify(g)
        s, t, h = gf_gcdex(F, G, mod, dom)
        return per(s), per(t), per(h)

    def gcd(f, g):
        r"""
        Returns polynomial GCD of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(8), ZZ(7)], 11, ZZ).gcd(
        ... GFP([ZZ(1), ZZ(7), ZZ(1), ZZ(7)], 11, ZZ)) == \
        ... GFP([ZZ(1), ZZ(7)], 11, ZZ)
        True
        """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_gcd(F, G, mod, dom))

    def lcm(f, g):
        r"""
        Returns polynomial LCM of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(8), ZZ(7)], 11, ZZ).lcm(
        ... GFP([ZZ(1), ZZ(7), ZZ(1), ZZ(7)], 11, ZZ)) == \
        ... GFP([ZZ(1), ZZ(8), ZZ(8), ZZ(8), ZZ(7)], 11, ZZ)
        True
        """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_lcm(F, G, mod, dom))

    def cofactors(f, g):
        r"""
        Returns polynomial GCD and cofactors of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(8), ZZ(7)], 11, ZZ).cofactors(
        ... GFP([ZZ(1), ZZ(7), ZZ(1), ZZ(7)], 11, ZZ)) == \
        ... (GFP([ZZ(1), ZZ(7)], 11, ZZ), GFP([ZZ(1), ZZ(1)], 11, ZZ),
        ...  GFP([ZZ(1), ZZ(0), ZZ(1)], 11, ZZ))
        True
        """
        mod, dom, per, F, G = f.unify(g)
        h, s, t = gf_cofactors(F, G, mod, dom)
        return per(h), per(s), per(t)

    def monic(f):
        r"""
        Divides all coefficients by `LC(f)`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).monic() == \
        ... GFP([ZZ(1), ZZ(4), ZZ(3)], 5, ZZ)
        True
        """
        return f.per(gf_monic(f.rep, f.mod, f.dom)[1])

    def diff(f, m=1):
        r"""
        Computes partial derivative of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).diff() == \
        ... GFP([ZZ(1), ZZ(2)], 5, ZZ)
        True
        """
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
        """
        Evaluates `f` at the given point `a`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).eval(2)
        0
        """
        return gf_eval(f.rep, f.dom.convert(a), f.mod, f.dom)

    def compose(f, g):
        r"""
        Computes functional composition of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).compose(
        ... GFP([ZZ(2), ZZ(2), ZZ(2)], 5, ZZ)) == \
        ... GFP([ZZ(2), ZZ(4), ZZ(0), ZZ(3), ZZ(0)], 5, ZZ)
        True
        """
        mod, dom, per, F, G = f.unify(g)
        return per(gf_compose(F, G, mod, dom))

    def sqf_part(f):
        r"""
        Computes square-free part of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(1), ZZ(3), ZZ(0), ZZ(1), ZZ(0), ZZ(2), ZZ(2), ZZ(1)],
        ... 5, ZZ).sqf_part() == \
        ... GFP([ZZ(1), ZZ(4), ZZ(3)], 5, ZZ)
        True
        """
        return f.per(gf_sqf_part(f.rep, f.mod, f.dom))

    def sqf_list(f):
        r"""
        Returns a list of square-free factors of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP({11: ZZ(1), 0: ZZ(1)}, 11, ZZ).sqf_list() == \
        ... (1, [(GFP([ZZ(1), ZZ(1)], 11, ZZ), 11)])
        True
        """
        coeff, factors = gf_sqf_list(f.rep, f.mod, f.dom)
        return coeff, [ (f.per(g), k) for g, k in factors ]

    def factor_list(f):
        r"""
        Returns a list of irreducible factors of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(5), ZZ(2), ZZ(7), ZZ(2)], 11, ZZ).factor_list() == \
        ... (5, [(GFP([ZZ(1), ZZ(2)], 11, ZZ), 1),
        ...  (GFP([ZZ(1), ZZ(8)], 11, ZZ), 2)])
        True
        """
        coeff, factors = gf_factor(f.rep, f.mod, f.dom)
        return coeff, [ (f.per(g), k) for g, k in factors ]

    @property
    def is_zero(f):
        """
        Returns `True` if `f` is a zero polynomial.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([], 5, ZZ).is_zero
        True
        >>> GFP([ZZ(1)], 5, ZZ).is_zero
        False
        """
        return not f

    @property
    def is_one(f):
        """
        Returns `True` if `f` is a unit polynomial.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([], 5, ZZ).is_one
        False
        >>> GFP([ZZ(1)], 5, ZZ).is_one
        True
        """
        return f.rep == [f.dom.one]

    @property
    def is_ground(f):
        """
        Returns `True` if `f` is an element of the ground domain.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(0)], 5, ZZ).is_ground
        False
        >>> GFP([ZZ(1)], 5, ZZ).is_ground
        True
        """
        return len(f.rep) <= 1

    @property
    def is_sqf(f):
        """
        Returns `True` if `f` is a square-free polynomial.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).is_sqf
        True
        >>> GFP([ZZ(2), ZZ(4), ZZ(4), ZZ(2), ZZ(2), ZZ(1), ZZ(4)], 5, ZZ).is_sqf
        False
        """
        return gf_sqf_p(f.rep, f.mod, f.dom)

    @property
    def is_monic(f):
        """
        Returns `True` if the leading coefficient of `f` is one.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).is_monic
        False
        >>> GFP([ZZ(1), ZZ(4), ZZ(3)], 5, ZZ).is_monic
        True
        """
        return f.dom.is_one(gf_LC(f.rep, f.dom))

    @property
    def is_irreducible(f):
        """
        Returns `True` if `f` has no factors over its domain.

        Example
        =======
        >>> from sympy.polys.polyclasses import GFP
        >>> from sympy.polys.algebratools import ZZ
        >>> GFP([ZZ(1), ZZ(4), ZZ(2), ZZ(2), ZZ(3), ZZ(2), ZZ(4), ZZ(1), ZZ(4),
        ... ZZ(0), ZZ(4)], 5, ZZ).is_irreducible
        True
        >>> GFP([ZZ(3), ZZ(2), ZZ(4)], 5, ZZ).is_irreducible
        False
        """
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
        """
        Unify representations of two multivariate polynomials.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ, QQ

        >> DMP([ZZ(1), ZZ(2)], ZZ).unify(DMP([[QQ(1)], [QQ(2)]], QQ))
        (0, ZZ, <bound method DMP.per of DMP([1, 2], ZZ)>, [1, 2], [[1], [2]])
        """
        # XXX: ???
        return f.lev, f.dom, f.per, f.rep, g.rep

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
        r"""
        Create a DMP out of the given representation.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> a = DMP([[ZZ(1)], [ZZ(2), ZZ(0)]], ZZ)
        >>> a.per([[ZZ(1), ZZ(0)], [ZZ(1)]]) == \
        ... DMP([[ZZ(1), ZZ(0)], [ZZ(1)]], ZZ)
        True
        >>> a.per([ZZ(1), ZZ(0)], ZZ, kill=True) == \
        ... DMP([ZZ(1), ZZ(0)], ZZ)
        True
        """
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
        """
        Returns a multivariate zero with level `lev` and domain `dom`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP.zero(3, ZZ)
        DMP([[[[]]]], ZZ)
        """
        return DMP(0, dom, lev)

    @classmethod
    def one(cls, lev, dom):
        r"""
        Returns a multivariate one with level `lev` and domain `dom`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP.one(3, ZZ) == \
        ... DMP([[[[ZZ(1)]]]], ZZ)
        True
        """
        return DMP(1, dom, lev)

    @classmethod
    def from_list(cls, rep, lev, dom):
        r"""
        Create an instance of `cls` given a list of native coefficients.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP.from_list([[1], [2, 3]], 1, ZZ) == \
        ... DMP([[ZZ(1)], [ZZ(2), ZZ(3)]], ZZ)
        True
        """
        return cls(dmp_convert(rep, lev, None, dom), dom, lev)

    @classmethod
    def from_sympy_list(cls, rep, lev, dom):
        r"""
        Create an instance of `cls` given a list of SymPy coefficients.

        Example
        =======
        >>> from sympy import S
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP.from_sympy_list([[S(1)], [S(2), S(3)]], 1, ZZ) == \
        ... DMP([[ZZ(1)], [ZZ(2), ZZ(3)]], ZZ)
        True
        """
        return cls(dmp_from_sympy(rep, lev, dom), dom, lev)

    def to_dict(f):
        """
        Convert `f` to a dict representation with native coefficients.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(1)], [ZZ(2), ZZ(0)]], ZZ).to_dict()
        {(0, 1): 2, (1, 0): 1}
        """
        return dmp_to_dict(f.rep, f.lev)

    def to_sympy_dict(f):
        """
        Convert `f` to a dict representation with SymPy coefficients.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(1)], [ZZ(2), ZZ(0)]], ZZ).to_sympy_dict()
        {(0, 1): 2, (1, 0): 1}
        """
        rep = dmp_to_dict(f.rep, f.lev)

        for k, v in rep.iteritems():
            rep[k] = f.dom.to_sympy(v)

        return rep

    def to_ring(f):
        r"""
        Make the ground domain a ring.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ, ZZ
        >>> DMP([[QQ(1)], [QQ(2), QQ(0)]], QQ).to_ring() == \
        ... DMP([[ZZ(1)], [ZZ(2), ZZ(0)]], ZZ)
        True
        """
        return f.convert(f.dom.get_ring())

    def to_field(f):
        r"""
        Make the ground domain a field.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ, QQ
        >>> DMP([[ZZ(1)], [ZZ(2), ZZ(0)]], ZZ).to_field() == \
        ... DMP([[QQ(1)], [QQ(2), QQ(0)]], QQ)
        True
        """
        return f.convert(f.dom.get_field())

    def to_exact(f):
        r"""
        Make the ground domain exact.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import RR, QQ
        >>> DMP([RR(1.5), RR(1.0)], RR).to_exact() == \
        ... DMP([QQ(3, 2), QQ(1)], QQ)
        True
        """
        return f.convert(f.dom.get_exact())

    def convert(f, dom):
        r"""
        Convert the ground domain of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ, QQ
        >>> DMP([ZZ(1), ZZ(2)], ZZ).convert(QQ) == \
        ... DMP([QQ(1), QQ(2)], QQ)
        True
        """
        if f.dom == dom:
            return f
        else:
            return DMP(dmp_convert(f.rep, f.lev, f.dom, dom), dom, f.lev)

    def coeffs(f, order=None):
        """
        Returns all non-zero coefficients from `f` in lex order.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(1)], [ZZ(2), ZZ(0), ZZ(-1)]], ZZ).coeffs()
        [1, 2, -1]
        """
        return [ c for _, c in dmp_list_terms(f.rep, f.lev, f.dom, order=order) ]

    def monoms(f, order=None):
        """
        Returns all non-zero monomials from `f` in lex order.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(1)], [ZZ(2), ZZ(0), ZZ(-1)]], ZZ).monoms()
        [(1, 0), (0, 2), (0, 0)]
        """
        return [ m for m, _ in dmp_list_terms(f.rep, f.lev, f.dom, order=order) ]

    def terms(f, order=None):
        """
        Returns all non-zero terms from `f` in lex order.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(1)], [ZZ(2), ZZ(0), ZZ(-1)]], ZZ).terms()
        [((1, 0), 1), ((0, 2), 2), ((0, 0), -1)]
        """
        return dmp_list_terms(f.rep, f.lev, f.dom, order=order)

    def all_coeffs(f):
        """
        Returns all coefficients from `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(2), ZZ(0), ZZ(-1)], ZZ).all_coeffs()
        [2, 0, -1]
        """
        if not f.lev:
            if not f:
                return [f.dom.zero]
            else:
                return [ c for c in f.rep ]
        else:
            raise PolynomialError('multivariate polynomials not supported')

    def all_monoms(f):
        """
        Returns all monomials from `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(2), ZZ(0), ZZ(-1)], ZZ).all_monoms()
        [(2,), (1,), (0,)]
        """
        if not f.lev:
            n = gf_degree(f.rep)

            if n < 0:
                return [(0,)]
            else:
                return [ (n-i,) for i, c in enumerate(f.rep) ]
        else:
            raise PolynomialError('multivariate polynomials not supported')

    def all_terms(f):
        """
        Returns all terms from a `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(2), ZZ(0), ZZ(-1)], ZZ).all_terms()
        [((2,), 2), ((1,), 0), ((0,), -1)]
        """
        if not f.lev:
            n = gf_degree(f.rep)

            if n < 0:
                return [((0,), f.dom.zero)]
            else:
                return [ ((n-i,), c) for i, c in enumerate(f.rep) ]
        else:
            raise PolynomialError('multivariate polynomials not supported')

    def lift(f):
        r"""
        Convert algebraic coefficients to rationals.

        Example
        =======
        >>> from sympy import I
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> K = QQ.algebraic_field(I)
        >>> DMP([K(1), K([QQ(1), QQ(0)]), K([QQ(2), QQ(0)])], K).lift() == \
        ... DMP([QQ(1), QQ(0), QQ(2), QQ(0), QQ(9), QQ(0), QQ(-8), QQ(0),
        ... QQ(16)], QQ)
        True
        """
        return f.per(dmp_lift(f.rep, f.lev, f.dom), dom=f.dom.dom)

    def deflate(f):
        r"""
        Reduce degree of `f` by mapping `x_i**m` to `y_i`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(0), ZZ(1), ZZ(0), ZZ(0), ZZ(1)],
        ... ZZ).deflate() == \
        ... ((3,), DMP([ZZ(1), ZZ(1), ZZ(1)], ZZ))
        True
        >>> DMP([[ZZ(1), ZZ(0), ZZ(0), ZZ(2)], [],
        ... [ZZ(3), ZZ(0), ZZ(0), ZZ(4)]], ZZ).deflate() == \
        ... ((2, 3), DMP([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], ZZ))
        True
        """
        J, F = dmp_deflate(f.rep, f.lev, f.dom)
        return J, f.per(F)

    def terms_gcd(f):
        r"""
        Remove GCD of terms from the polynomial `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(1), ZZ(0)], [ZZ(1), ZZ(0), ZZ(0)], [], []],
        ... ZZ).terms_gcd() == \
        ... ((2, 1), DMP([[ZZ(1)], [ZZ(1), ZZ(0)]], ZZ))
        True
        """
        J, F = dmp_terms_gcd(f.rep, f.lev, f.dom)
        return J, f.per(F)

    def mul_ground(f, c):
        r"""
        Multiply `f` by a an element of the ground domain.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(2)], [ZZ(2), ZZ(0)]], ZZ).mul_ground(3) == \
        ... DMP([[ZZ(6)], [ZZ(6), ZZ(0)]], ZZ)
        True
        """
        return f.per(dmp_mul_ground(f.rep, f.dom.convert(c), f.lev, f.dom))

    def quo_ground(f, c):
        r"""
        Quotient of `f` by a an element of the ground domain.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> DMP([[QQ(1), QQ(0)], [QQ(2)], []], QQ).quo_ground(QQ(2)) == \
        ... DMP([[QQ(1, 2), QQ(0)], [QQ(1)], []], QQ)
        True
        """
        return f.per(dmp_quo_ground(f.rep, f.dom.convert(c), f.lev, f.dom))

    def exquo_ground(f, c):
        r"""
        Exact quotient of `f` by a an element of the ground domain.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ, QQ
        >>> DMP([ZZ(3), ZZ(0), ZZ(2)], ZZ).exquo_ground(ZZ(2)) == \
        ... DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ)
        True
        >>> DMP([QQ(3), QQ(0), QQ(2)], QQ).exquo_ground(QQ(2)) == \
        ... DMP([QQ(3, 2), QQ(0), QQ(1)], QQ)
        True
        """
        return f.per(dmp_exquo_ground(f.rep, f.dom.convert(c), f.lev, f.dom))

    def abs(f):
        r"""
        Make all coefficients in `f` positive.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(1), ZZ(0)], [ZZ(-1)], []], ZZ).abs() == \
        ... DMP([[ZZ(1), ZZ(0)], [ZZ(1)], []], ZZ)
        True
        """
        return f.per(dmp_abs(f.rep, f.lev, f.dom))

    def neg(f):
        r"""
        Negate all cefficients in `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(1), ZZ(0)], [ZZ(-1)], []], ZZ).neg() == \
        ... DMP([[ZZ(-1), ZZ(0)], [ZZ(1)], []], ZZ)
        True
        """
        return f.per(dmp_neg(f.rep, f.lev, f.dom))

    def add(f, g):
        r"""
        Add two multivariate polynomials `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ).add(DMP([ZZ(1), ZZ(-2)], ZZ)) == \
        ... DMP([ZZ(1), ZZ(1), ZZ(-3)], ZZ)
        True
        """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_add(F, G, lev, dom))

    def sub(f, g):
        r"""
        Subtract two multivariate polynomials `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ).sub(DMP([ZZ(1), ZZ(-2)], ZZ)) == \
        ... DMP([ZZ(1), ZZ(-1), ZZ(1)], ZZ)
        True
        """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_sub(F, G, lev, dom))

    def mul(f, g):
        r"""
        Multiply two multivariate polynomials `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(2)], ZZ).mul(DMP([ZZ(1), ZZ(-2)], ZZ)) == \
        ... DMP([ZZ(1), ZZ(0), ZZ(-4)], ZZ)
        True
        """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_mul(F, G, lev, dom))

    def sqr(f):
        r"""
        Square a multivariate polynomial `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ).sqr() == \
        ... DMP([ZZ(1), ZZ(0), ZZ(2), ZZ(0), ZZ(1)], ZZ)
        True
        """
        return f.per(dmp_sqr(f.rep, f.lev, f.dom))

    def pow(f, n):
        r"""
        Raise `f` to a non-negative power `n`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(-2)], ZZ).pow(3) == \
        ... DMP([ZZ(1), ZZ(-6), ZZ(12), ZZ(-8)], ZZ)
        True
        """
        if isinstance(n, int):
            return f.per(dmp_pow(f.rep, n, f.lev, f.dom))
        else:
            raise TypeError("`int` expected, got %s" % type(n))

    def pdiv(f, g):
        r"""
        Polynomial pseudo-division of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ).pdiv(DMP([ZZ(2), ZZ(-4)], ZZ)) == \
        ... (DMP([ZZ(2), ZZ(4)], ZZ), DMP([ZZ(20)], ZZ))
        True
        """
        lev, dom, per, F, G = f.unify(g)
        q, r = dmp_pdiv(F, G, lev, dom)
        return per(q), per(r)

    def prem(f, g):
        r"""
        Polynomial pseudo-remainder of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ).prem(DMP([ZZ(2), ZZ(-4)], ZZ)) == \
        ... DMP([ZZ(20)], ZZ)
        True
        """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_prem(F, G, lev, dom))

    def pquo(f, g):
        r"""
        Polynomial pseudo-quotient of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ

        >> DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ).pquo(DMP([ZZ(2), ZZ(-4)], ZZ))
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: [2, -4] does not divide [1, 0, 1]
        >>> DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ).pquo(DMP([ZZ(2), ZZ(-2)], ZZ)) == \
        ... DMP([ZZ(2), ZZ(2)], ZZ)
        True
        """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_pquo(F, G, lev, dom))

    def pexquo(f, g):
        r"""
        Polynomial exact pseudo-quotient of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ).pexquo(DMP([ZZ(2), ZZ(-4)], ZZ)) == \
        ... DMP([ZZ(2), ZZ(4)], ZZ)
        True
        >>> DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ).pexquo(DMP([ZZ(2), ZZ(-2)], ZZ)) == \
        ... DMP([ZZ(2), ZZ(2)], ZZ)
        True
        """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_pexquo(F, G, lev, dom))

    def div(f, g):
        r"""
        Polynomial division with remainder of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> DMP([QQ(1), QQ(0), QQ(1)], QQ).div(DMP([QQ(2), QQ(-4)], QQ)) == \
        ... (DMP([QQ(1, 2), QQ(1)], QQ), DMP([QQ(5)], QQ))
        True
        """
        lev, dom, per, F, G = f.unify(g)
        q, r = dmp_div(F, G, lev, dom)
        return per(q), per(r)

    def rem(f, g):
        r"""
        Computes polynomial remainder of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> DMP([QQ(1), QQ(0), QQ(1)], QQ).rem(DMP([QQ(2), QQ(-4)], QQ)) == \
        ... DMP([QQ(5)], QQ)
        True
        """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_rem(F, G, lev, dom))

    def quo(f, g):
        r"""
        Computes polynomial quotient of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ

        >> DMP([QQ(1), QQ(0), QQ(1)], QQ).quo(DMP([QQ(2), QQ(-4)], QQ))
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: [1, -1] does not divide [1, 0, 1]
        >>> DMP([QQ(1), QQ(0), QQ(-1)], QQ).quo(DMP([QQ(1), QQ(-1)], QQ)) == \
        ... DMP([QQ(1), QQ(1)], QQ)
        True
        """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_quo(F, G, lev, dom))

    def exquo(f, g):
        r"""

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> DMP([QQ(1), QQ(0), QQ(1)], QQ).exquo(DMP([QQ(2), QQ(-4)], QQ)) == \
        ... DMP([QQ(1, 2), QQ(1)], QQ)
        True
        >>> DMP([QQ(1), QQ(0), QQ(-1)], QQ).exquo(DMP([QQ(1), QQ(-1)], QQ)) == \
        ... DMP([QQ(1), QQ(1)], QQ)
        True
        """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_exquo(F, G, lev, dom))

    def degree(f, j=0):
        """
        Returns the leading degree of `f` in `x_j`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(2)], [ZZ(1), ZZ(2), ZZ(3)]], ZZ).degree()
        1
        >>> DMP([[ZZ(2)], [ZZ(1), ZZ(2), ZZ(3)]], ZZ).degree(1)
        2
        """
        if isinstance(j, int):
            return dmp_degree_in(f.rep, j, f.lev)
        else:
            raise TypeError("`int` expected, got %s" % type(j))

    def degree_list(f):
        """
        Returns a list of degrees of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(2)], [ZZ(1), ZZ(2), ZZ(3)]], ZZ).degree_list()
        (1, 2)
        """
        return dmp_degree_list(f.rep, f.lev)

    def total_degree(f):
        """
        Returns the total degree of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(2)], [ZZ(1), ZZ(2), ZZ(3)]], ZZ).total_degree()
        3
        """
        return sum(dmp_degree_list(f.rep, f.lev))

    def LC(f):
        """
        Returns the leading coefficent of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(1)], [ZZ(2), ZZ(3)]], ZZ).LC()
        1
        """
        return dmp_ground_LC(f.rep, f.lev, f.dom)

    def TC(f):
        """
        Returns the trailing coefficent of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(1)], [ZZ(2), ZZ(3)]], ZZ).TC()
        3
        """
        return dmp_ground_TC(f.rep, f.lev, f.dom)

    def nth(f, *N):
        """
        Returns the `n`-th coefficient of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(1)], [ZZ(2), ZZ(3)]], ZZ).nth(0, 1)
        2
        """
        if all(isinstance(n, (int, long)) for n in N):
            return dmp_ground_nth(f.rep, N, f.lev, f.dom)
        else:
            raise TypeError("a sequence of integers expected")

    def max_norm(f):
        """
        Returns maximum norm of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(-1), ZZ(2), ZZ(3)], ZZ).max_norm()
        3
        """
        return dmp_max_norm(f.rep, f.lev, f.dom)

    def l1_norm(f):
        """
        Returns l1 norm of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(2), ZZ(-3), ZZ(0), ZZ(1)], ZZ).l1_norm()
        6
        """
        return dmp_l1_norm(f.rep, f.lev, f.dom)

    def clear_denoms(f):
        r"""
        Clear denominators, but keep the ground domain.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ, ZZ
        >>> DMP([QQ(1, 2), QQ(1, 3)], QQ).clear_denoms() == \
        ... (ZZ(6), DMP([QQ(3), QQ(2)], QQ))
        True
        """
        coeff, F = dmp_clear_denoms(f.rep, f.lev, f.dom)
        return coeff, f.per(F)

    def integrate(f, m=1, j=0):
        r"""
        Computes `m`-th order indefinite integral of `f` in `x_j`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> DMP([QQ(1), QQ(2), QQ(0)], QQ).integrate() == \
        ... DMP([QQ(1, 3), QQ(1), QQ(0), QQ(0)], QQ)
        True
        >>> DMP([QQ(1), QQ(2), QQ(0)], QQ).integrate(m=2) == \
        ... DMP([QQ(1, 12), QQ(1, 3), QQ(0), QQ(0), QQ(0)], QQ)
        True
        >>> DMP([[QQ(1)], [QQ(2), QQ(0)]], QQ).integrate(j=0) == \
        ... DMP([[QQ(1, 2)], [QQ(2), QQ(0)], []], QQ)
        True
        """
        if not isinstance(m, int):
            raise TypeError("`int` expected, got %s" % type(m))

        if not isinstance(j, int):
            raise TypeError("`int` expected, got %s" % type(j))

        return f.per(dmp_integrate_in(f.rep, m, j, f.lev, f.dom))

    def diff(f, m=1, j=0):
        r"""
        Computes `m`-th order derivative of `f` in `x_j`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(2), ZZ(3), ZZ(4)], ZZ).diff() == \
        ... DMP([ZZ(3), ZZ(4), ZZ(3)], ZZ)
        True
        >>> DMP([ZZ(1), ZZ(2), ZZ(3), ZZ(4)], ZZ).diff(m=2) == \
        ... DMP([ZZ(6), ZZ(4)], ZZ)
        True
        >>> DMP([[1, 2, 3], [2, 3, 1]], ZZ).diff(j=1) == \
        ... DMP([[ZZ(2), ZZ(2)], [ZZ(4), ZZ(3)]], ZZ)
        True
        """
        if not isinstance(m, int):
            raise TypeError("`int` expected, got %s" % type(m))

        if not isinstance(j, int):
            raise TypeError("`int` expected, got %s" % type(j))

        return f.per(dmp_diff_in(f.rep, m, j, f.lev, f.dom))

    def eval(f, a, j=0):
        r"""
        Evaluates `f` at the given point `a` in `x_j`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(2), ZZ(3)], [ZZ(1), ZZ(2)]], ZZ).eval(2) == \
        ... DMP([ZZ(5), ZZ(8)], ZZ)
        True
        >>> DMP([[ZZ(2), ZZ(3)], [ZZ(1), ZZ(2)]], ZZ).eval(2, j=1) == \
        ... DMP([ZZ(7), ZZ(4)], ZZ)
        True
        """
        if not isinstance(j, int):
            raise TypeError("`int` expected, got %s" % type(j))

        return f.per(dmp_eval_in(f.rep,
            f.dom.convert(a), j, f.lev, f.dom), kill=True)

    def half_gcdex(f, g):
        r"""
        Half extended Euclidean algorithm, if univariate.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> DMP([QQ(1), QQ(1), QQ(0)], QQ).half_gcdex(
        ... DMP([QQ(1), QQ(3), QQ(2)], QQ)) == \
        ... (DMP([QQ(-1, 2)], QQ), DMP([QQ(1), QQ(1)], QQ))
        True
        """
        lev, dom, per, F, G = f.unify(g)

        if not lev:
            s, h = dup_half_gcdex(F, G, dom)
            return per(s), per(h)
        else:
            raise ValueError('univariate polynomial expected')

    def gcdex(f, g):
        r"""
        Extended Euclidean algorithm, if univariate.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> DMP([QQ(1), QQ(1), QQ(0)], QQ).gcdex(
        ... DMP([QQ(1), QQ(3), QQ(2)], QQ)) == \
        ... (DMP([QQ(-1, 2)], QQ), DMP([QQ(1, 2)], QQ), DMP([QQ(1), QQ(1)], QQ))
        True
        """
        lev, dom, per, F, G = f.unify(g)

        if not lev:
            s, t, h = dup_gcdex(F, G, dom)
            return per(s), per(t), per(h)
        else:
            raise ValueError('univariate polynomial expected')

    def invert(f, g):
        r"""
        Invert `f` modulo `g`, if possible.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> DMP([QQ(1), QQ(0), QQ(-1)], QQ).invert(DMP([QQ(2), QQ(-1)], QQ)) == \
        ... DMP([QQ(-4, 3)], QQ)
        True
        >>> DMP([QQ(1), QQ(0), QQ(-1)], QQ).invert(DMP([QQ(1), QQ(-1)], QQ))
        Traceback (most recent call last):
        ...
        NotInvertible: zero divisor
        """
        lev, dom, per, F, G = f.unify(g)

        if not lev:
            return per(dup_invert(F, G, dom))
        else:
            raise ValueError('univariate polynomial expected')

    def subresultants(f, g):
        r"""
        Computes subresultant PRS sequence of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ).subresultants(
        ... DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ)) == \
        ... [DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ), DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ),
        ...  DMP([ZZ(-2)], ZZ)]
        True
        """
        lev, dom, per, F, G = f.unify(g)
        R = dmp_subresultants(F, G, lev, dom)
        return map(per, R)

    def resultant(f, g):
        r"""
        Computes resultant of `f` and `g` via PRS.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ).resultant(
        ... DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ))
        4
        """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_resultant(F, G, lev, dom), kill=True)

    def discriminant(f):
        """
        Computes discriminant of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(2), ZZ(3)], ZZ).discriminant()
        -8
        """
        return f.per(dmp_discriminant(f.rep, f.lev, f.dom), kill=True)

    def cofactors(f, g):
        r"""
        Returns GCD of `f` and `g` and their cofactors.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ).cofactors(
        ... DMP([ZZ(1), ZZ(-3), ZZ(2)], ZZ)) == \
        ... (DMP([ZZ(1), ZZ(-1)], ZZ), DMP([ZZ(1), ZZ(1)], ZZ),
        ...  DMP([ZZ(1), ZZ(-2)], ZZ))
        True
        """
        lev, dom, per, F, G = f.unify(g)
        h, cff, cfg = dmp_inner_gcd(F, G, lev, dom)
        return per(h), per(cff), per(cfg)

    def gcd(f, g):
        r"""
        Returns polynomial GCD of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ).gcd(
        ... DMP([ZZ(1), ZZ(-3), ZZ(2)], ZZ)) == \
        ... DMP([ZZ(1), ZZ(-1)], ZZ)
        True
        """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_gcd(F, G, lev, dom))

    def lcm(f, g):
        r"""
        Returns polynomial LCM of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ).lcm(
        ... DMP([ZZ(1), ZZ(-3), ZZ(2)], ZZ)) == \
        ... DMP([ZZ(1), ZZ(-2), ZZ(-1), ZZ(2)], ZZ)
        True
        """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_lcm(F, G, lev, dom))

    def trunc(f, p):
        r"""
        Reduce `f` modulo a constant `p`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(2), ZZ(3), ZZ(5), ZZ(7)], ZZ).trunc(3) == \
        ... DMP([ZZ(-1), ZZ(0), ZZ(-1), ZZ(1)], ZZ)
        True
        """
        return f.per(dmp_ground_trunc(f.rep, f.dom.convert(p), f.lev, f.dom))

    def monic(f):
        r"""
        Divides all coefficients by `LC(f)`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ, ZZ
        >>> DMP([ZZ(3), ZZ(6), ZZ(9)], ZZ).monic() == \
        ... DMP([ZZ(1), ZZ(2), ZZ(3)], ZZ)
        True
        >>> DMP([QQ(3), QQ(4), QQ(2)], QQ).monic() == \
        ... DMP([QQ(1), QQ(4, 3), QQ(2, 3)], QQ)
        True
        """
        return f.per(dmp_ground_monic(f.rep, f.lev, f.dom))

    def content(f):
        """
        Returns GCD of polynomial coefficients.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(2), ZZ(6)], [ZZ(4), ZZ(12)]], ZZ).content()
        2
        """
        return dmp_ground_content(f.rep, f.lev, f.dom)

    def primitive(f):
        r"""
        Returns content and a primitive form of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(2), ZZ(6)], [ZZ(4), ZZ(12)]], ZZ).primitive() == \
        ... (2, DMP([[ZZ(1), ZZ(3)], [ZZ(2), ZZ(6)]], ZZ))
        True
        """
        cont, F = dmp_ground_primitive(f.rep, f.lev, f.dom)
        return cont, f.per(F)

    def compose(f, g):
        r"""
        Computes functional composition of `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(1), ZZ(0)], ZZ).compose(DMP([ZZ(1), ZZ(-1)], ZZ)) == \
        ... DMP([ZZ(1), ZZ(-1), ZZ(0)], ZZ)
        True
        """
        lev, dom, per, F, G = f.unify(g)
        return per(dmp_compose(F, G, lev, dom))

    def decompose(f):
        r"""
        Computes functional decomposition of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(-2), ZZ(1), ZZ(0), ZZ(0)], ZZ).decompose() == \
        ... [DMP([ZZ(1), ZZ(0), ZZ(0)], ZZ), DMP([ZZ(1), ZZ(-1), ZZ(0)], ZZ)]
        True
        """
        if not f.lev:
            return map(f.per, dup_decompose(f.rep, f.dom))
        else:
            raise ValueError('univariate polynomial expected')

    def sturm(f):
        r"""
        Computes the Sturm sequence of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> DMP([QQ(1), QQ(-2), QQ(1), QQ(-3)], QQ).sturm() == \
        ... [DMP([QQ(1), QQ(-2), QQ(1), QQ(-3)], QQ),
        ...  DMP([QQ(3), QQ(-4), QQ(1)], QQ),
        ...  DMP([QQ(2, 9), QQ(25, 9)], QQ), DMP([QQ(-2079, 4)], QQ)]
        True
        """
        if not f.lev:
            return map(f.per, dup_sturm(f.rep, f.dom))
        else:
            raise ValueError('univariate polynomial expected')

    def sqf_norm(f):
        r"""
        Computes square-free norm of `f`.

        Example
        =======
        >>> from sympy import sqrt
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import QQ
        >>> K = QQ.algebraic_field(sqrt(3))
        >>> DMP([K(1), K(0), K(-2)], K).sqf_norm() == \
        ... (1, DMP([K(1), K([QQ(-2), QQ(0)]), K(1)], K),
        ...  DMP([QQ(1), QQ(0), QQ(-10), QQ(0), QQ(1)], QQ))
        True
        """
        s, g, r = dmp_sqf_norm(f.rep, f.lev, f.dom)
        return s, f.per(g), f.per(r, dom=f.dom.dom)

    def sqf_part(f):
        r"""
        Computes square-free part of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(-3), ZZ(-2)], ZZ).sqf_part() == \
        ... DMP([ZZ(1), ZZ(-1), ZZ(-2)], ZZ)
        True
        """
        return f.per(dmp_sqf_part(f.rep, f.lev, f.dom))

    def sqf_list(f, all=False):
        r"""
        Returns a list of square-free factors of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(2), ZZ(16), ZZ(50), ZZ(76), ZZ(56), ZZ(16)],
        ... ZZ).sqf_list() == \
        ... (2, [(DMP([ZZ(1), ZZ(1)], ZZ), 2), (DMP([ZZ(1), ZZ(2)], ZZ), 3)])
        True
        >>> DMP([ZZ(2), ZZ(16), ZZ(50), ZZ(76), ZZ(56), ZZ(16)],
        ... ZZ).sqf_list(all=True) == \
        ... (2, [(DMP([ZZ(1)], ZZ), 1), (DMP([ZZ(1), ZZ(1)], ZZ), 2),
        ... (DMP([ZZ(1), ZZ(2)], ZZ), 3)])
        True
        """
        coeff, factors = dmp_sqf_list(f.rep, f.lev, f.dom, all)
        return coeff, [ (f.per(g), k) for g, k in factors ]

    def sqf_list_include(f, all=False):
        r"""
        Returns a list of square-free factors of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(2), ZZ(16), ZZ(50), ZZ(76), ZZ(56), ZZ(16)],
        ... ZZ).sqf_list_include() == \
        ... [(DMP([ZZ(2), ZZ(2)], ZZ), 2), (DMP([ZZ(1), ZZ(2)], ZZ), 3)]
        True
        >>> DMP([ZZ(2), ZZ(16), ZZ(50), ZZ(76), ZZ(56), ZZ(16)],
        ... ZZ).sqf_list_include(all=True) == \
        ... [(DMP([ZZ(2)], ZZ), 1), (DMP([ZZ(1), ZZ(1)], ZZ), 2),
        ... (DMP([ZZ(1), ZZ(2)], ZZ), 3)]
        True
        """
        factors = dmp_sqf_list_include(f.rep, f.lev, f.dom, all)
        return [ (f.per(g), k) for g, k in factors ]

    def factor_list(f):
        r"""
        Returns a list of irreducible factors of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(2), ZZ(2), ZZ(2), ZZ(2), ZZ(0), ZZ(0)],
        ... ZZ).factor_list() == \
        ... (2, [(DMP([ZZ(1), ZZ(1)], ZZ), 1), (DMP([ZZ(1), ZZ(0)], ZZ), 2),
        ... (DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ), 1)])
        True
        """
        coeff, factors = dmp_factor_list(f.rep, f.lev, f.dom)
        return coeff, [ (f.per(g), k) for g, k in factors ]

    def factor_list_include(f):
        r"""
        Returns a list of irreducible factors of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(2), ZZ(2), ZZ(2), ZZ(2), ZZ(0), ZZ(0)],
        ... ZZ).factor_list_include() == \
        ... [(DMP([ZZ(2), ZZ(2)], ZZ), 1), (DMP([ZZ(1), ZZ(0)], ZZ), 2),
        ... (DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ), 1)]
        True
        """
        factors = dmp_factor_list_include(f.rep, f.lev, f.dom)
        return [ (f.per(g), k) for g, k in factors ]

    def intervals(f, all=False, eps=None, inf=None, sup=None, fast=False, sqf=False):
        r"""
        Compute isolating intervals for roots of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ, QQ
        >>> DMP([ZZ(1), ZZ(0), ZZ(-2)], ZZ).intervals() == \
        ... [((QQ(-2), QQ(-1)), 1), ((QQ(1), QQ(2)), 1)]
        True
        """
        if not f.lev:
            if not all:
                if not sqf:
                    return dup_isolate_real_roots(f.rep, f.dom, eps=eps, inf=inf, sup=sup, fast=fast)
                else:
                    return dup_isolate_real_roots_sqf(f.rep, f.dom, eps=eps, inf=inf, sup=sup, fast=fast)
            else:
                if not sqf:
                    return dup_isolate_all_roots(f.rep, f.dom, eps=eps, inf=inf, sup=sup, fast=fast)
                else:
                    return dup_isolate_all_roots_sqf(f.rep, f.dom, eps=eps, inf=inf, sup=sup, fast=fast)
        else:
            raise PolynomialError("can't isolate roots of a multivariate polynomial")

    def refine_root(f, s, t, eps=None, steps=None, fast=False):
        r"""
        Refine an isolating interval to the given precision.

        Example
        =======
        >>> from sympy import S
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ, QQ
        >>> DMP([ZZ(1), ZZ(0), ZZ(-2)], ZZ).refine_root(QQ(1), QQ(2),
        ... eps=1e-2) == \
        ... (QQ(24, 17), QQ(17, 12))
        True
        """
        if not f.lev:
            return dup_refine_real_root(f.rep, s, t, f.dom, eps=eps, steps=steps, fast=fast)
        else:
            raise PolynomialError("can't refine a root of a multivariate polynomial")

    def count_real_roots(f, inf=None, sup=None):
        """
        Return the number of real roots of ``f`` in ``[inf, sup]``.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(0), ZZ(0), ZZ(-4)], ZZ).count_real_roots()
        2
        """
        return dup_count_real_roots(f.rep, f.dom, inf=inf, sup=sup)

    def count_complex_roots(f, inf=None, sup=None):
        """
        Return the number of complex roots of ``f`` in ``[inf, sup]``.

        Example
        =======
        >>> from sympy import I
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(0), ZZ(0), ZZ(0), ZZ(-4)], ZZ).count_complex_roots()
        4
        """
        return dup_count_complex_roots(f.rep, f.dom, inf=inf, sup=sup)

    @property
    def is_zero(f):
        """
        Returns `True` if `f` is a zero polynomial.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[[[[[[]]]]]]], ZZ).is_zero
        True
        >>> DMP([[[[[[[ZZ(1)]]]]]]], ZZ).is_zero
        False
        """
        return dmp_zero_p(f.rep, f.lev)

    @property
    def is_one(f):
        """
        Returns `True` if `f` is a unit polynomial.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[[[[[[]]]]]]], ZZ).is_one
        False
        >>> DMP([[[[[[[ZZ(1)]]]]]]], ZZ).is_one
        True
        """
        return dmp_one_p(f.rep, f.lev, f.dom)

    @property
    def is_ground(f):
        """
        Returns `True` if `f` is an element of the ground domain.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[[ZZ(3)]]], ZZ).is_ground
        True
        >>> DMP([[ZZ(1)], []], ZZ).is_ground
        False
        """
        return dmp_ground_p(f.rep, None, f.lev)

    @property
    def is_sqf(f):
        """
        Returns `True` if `f` is a square-free polynomial.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(-2), ZZ(1)], ZZ).is_sqf
        False
        >>> DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ).is_sqf
        True
        """
        return dmp_sqf_p(f.rep, f.lev, f.dom)

    @property
    def is_monic(f):
        """
        Returns `True` if the leading coefficient of `f` is one.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(1), ZZ(2)], ZZ).is_monic
        True
        >>> DMP([ZZ(2), ZZ(2)], ZZ).is_monic
        False
        """
        return f.dom.is_one(dmp_ground_LC(f.rep, f.lev, f.dom))

    @property
    def is_primitive(f):
        """
        Returns `True` if GCD of coefficients of `f` is one.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([ZZ(2), ZZ(6), ZZ(12)], ZZ).is_primitive
        False
        >>> DMP([ZZ(1), ZZ(3), ZZ(6)], ZZ).is_primitive
        True
        """
        return f.dom.is_one(dmp_ground_content(f.rep, f.lev, f.dom))

    @property
    def is_linear(f):
        """
        Returns `True` if `f` is linear in all its variables.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(1), ZZ(0)], [ZZ(2)]], ZZ).is_linear
        True
        >>> DMP([[ZZ(1), ZZ(0)], [ZZ(2)], []], ZZ).is_linear
        False
        """
        return all([ sum(monom) <= 1 for monom in dmp_to_dict(f.rep, f.lev).keys() ])

    @property
    def is_homogeneous(f):
        """
        Returns `True` if `f` has zero trailing coefficient.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMP([[ZZ(1), ZZ(1)], [ZZ(1), ZZ(0)]], ZZ).is_homogeneous
        True
        >>> DMP([[ZZ(1), ZZ(1)], [ZZ(1), ZZ(1)]], ZZ).is_homogeneous
        False
        """
        # XXX: This is not the same as the definition of homogeneous at
        # http://en.wikipedia.org/wiki/Homogeneous_polynomial.
        # c.f. also the homogeneous_order() function in ode.py.
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
        r"""
        Returns a multivariate zero fraction with level `lev` and domain `dom`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF.zero(3, ZZ) == \
        ... DMF(([[[[]]]], [[[[ZZ(1)]]]]), ZZ)
        True
        """
        return DMF(0, dom, lev)

    @classmethod
    def one(cls, lev, dom):
        r"""
        Returns a multivariate zero fraction with level `lev` and domain `dom`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF.one(3, ZZ) == \
        ... DMF(([[[[ZZ(1)]]]], [[[[ZZ(1)]]]]), ZZ)
        True
        """
        return DMF(1, dom, lev)

    def numer(f):
        r"""
        Returns the numerator of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMF, DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF(([ZZ(1), ZZ(2)], [ZZ(1), ZZ(-2)]), ZZ).numer() == \
        ... DMP([ZZ(1), ZZ(2)], ZZ)
        True
        """
        return f.half_per(f.num)

    def denom(f):
        r"""
        Returns the denominator of `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMF, DMP
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF(([ZZ(1), ZZ(2)], [ZZ(1), ZZ(-2)]), ZZ).denom() == \
        ... DMP([ZZ(1), ZZ(-2)], ZZ)
        True
        """
        return f.half_per(f.den)

    def cancel(f):
        r"""
        Remove common factors from `f.num` and `f.den`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF(([ZZ(1), ZZ(0), ZZ(-1)], [ZZ(1), ZZ(1)]), ZZ).cancel() == \
        ... DMF(([ZZ(1), ZZ(-1)], [ZZ(1)]), ZZ)
        True
        """
        return f.per(f.num, f.den)

    def neg(f):
        r"""
        Negate all cefficients in `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF(([ZZ(1), ZZ(0), ZZ(-1)], [ZZ(1), ZZ(-2)]), ZZ).neg() == \
        ... DMF(([ZZ(-1), ZZ(0), ZZ(1)], [ZZ(1), ZZ(-2)]), ZZ)
        True
        """
        return f.per(dmp_neg(f.num, f.lev, f.dom), f.den, cancel=False)

    def add(f, g):
        r"""
        Add two multivariate fractions `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF(([ZZ(1), ZZ(1)], [ZZ(1), ZZ(-1)]), ZZ).add(
        ... DMF(([ZZ(1), ZZ(2)], [ZZ(1)]), ZZ)) == \
        ... DMF(([ZZ(1), ZZ(2), ZZ(-1)], [ZZ(1), ZZ(-1)]), ZZ)
        True
        """
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
        r"""
        Subtract two multivariate fractions `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF(([ZZ(1), ZZ(1)], [ZZ(1), ZZ(-1)]), ZZ).sub(
        ... DMF(([ZZ(1), ZZ(2)], [ZZ(1)]), ZZ)) == \
        ... DMF(([ZZ(-1), ZZ(0), ZZ(3)], [ZZ(1), ZZ(-1)]), ZZ)
        True
        """
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
        r"""
        Multiply two multivariate fractions `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF(([ZZ(1), ZZ(1)], [ZZ(1), ZZ(-1)]), ZZ).mul(
        ... DMF(([ZZ(1), ZZ(2)], [ZZ(2), ZZ(-1)]), ZZ)) == \
        ... DMF(([ZZ(1), ZZ(3), ZZ(2)], [ZZ(2), ZZ(-3), ZZ(1)]), ZZ)
        True
        """
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
        r"""
        Raise `f` to a non-negative power `n`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF(([ZZ(1), ZZ(-1)], [ZZ(1), ZZ(-2)]), ZZ).pow(2) == \
        ... DMF(([ZZ(1), ZZ(-2), ZZ(1)], [ZZ(1), ZZ(-4), ZZ(4)]), ZZ)
        True
        """
        if isinstance(n, int):
            return f.per(dmp_pow(f.num, n, f.lev, f.dom),
                         dmp_pow(f.den, n, f.lev, f.dom), cancel=False)
        else:
            raise TypeError("`int` expected, got %s" % type(n))

    def quo(f, g):
        r"""
        Computes quotient of fractions `f` and `g`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF(([ZZ(1), ZZ(-1)], [ZZ(1), ZZ(2)]), ZZ).quo(
        ... DMF(([ZZ(1), ZZ(2)], [ZZ(2), ZZ(-1)]), ZZ)) == \
        ... DMF(([ZZ(2), ZZ(-3), ZZ(1)], [ZZ(1), ZZ(4), ZZ(4)]), ZZ)
        True
        """
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
        r"""
        Computes inverse of a fraction `f`.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF(([ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]), ZZ).invert() == \
        ... DMF(([ZZ(3), ZZ(4)], [ZZ(1), ZZ(2)]), ZZ)
        True
        """
        return f.per(f.den, f.num, cancel=False)

    @property
    def is_zero(f):
        """
        Returns `True` if `f` is a zero fraction.

        Example
        =======
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF(([[]], [[ZZ(1)]]), ZZ).is_zero
        True
        >>> DMF(([[ZZ(1)]], [ZZ(1)]), ZZ).is_zero
        False
        """
        return dmp_zero_p(f.num, f.lev)

    @property
    def is_one(f):
        """
        Returns `True` if `f` is a unit fraction.
        Example
        =======
        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.algebratools import ZZ
        >>> DMF(([[]], [[ZZ(1)]]), ZZ).is_one
        False
        >>> DMF(([[ZZ(1)]], [ZZ(1)]), ZZ).is_one
        True
        """
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

        if isinstance(mod, DMP):
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

