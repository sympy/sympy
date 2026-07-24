"""Low-level sparse polynomial operations on raw dict representations.

These helpers deliberately do not depend on PolyRing or PolyElement.
They operate on dicts mapping exponent tuples to coefficients and use
only the provided domain and number of variables.
"""

from __future__ import annotations

from itertools import combinations
from typing import TYPE_CHECKING, Mapping, Sequence, TypeVar

from sympy.core.intfunc import igcd
from sympy.ntheory.multinomial import multinomial_coefficients
from sympy.polys.monomials import monom
from sympy.polys.orderings import lex, MonomialOrder

if TYPE_CHECKING:
    from sympy.polys.domains.domain import Domain, Er

_T = TypeVar("_T")

smp = dict[monom, _T]


def smp_add(f: smp[Er], g: smp[Er], domain: Domain[Er], n: int) -> smp[Er]:
    # Add two sparse polynomials and return a new dictionary.
    zero = domain.zero
    h = f.copy()
    for mon, coeff in g.items():
        coeff = h.get(mon, zero) + coeff
        if coeff:
            h[mon] = coeff
        elif mon in h:
            del h[mon]
    return h


def smp_mul(f: smp[Er], g: smp[Er], domain: Domain[Er], n: int) -> smp[Er]:
    # Multiply two sparse polynomials and return a new dictionary.
    zero = domain.zero
    h: smp[Er] = {}
    for mon1, coeff1 in f.items():
        for mon2, coeff2 in g.items():
            mon = tuple(mon1[i] + mon2[i] for i in range(n))
            coeff = h.get(mon, zero) + coeff1 * coeff2
            if coeff:
                h[mon] = coeff
            elif mon in h:
                del h[mon]
    return h


def smp_sub(f: smp[Er], g: smp[Er], domain: Domain[Er], n: int) -> smp[Er]:
    # Subtract one sparse polynomial from another.
    zero = domain.zero
    h = f.copy()
    for mon, coeff in g.items():
        coeff = h.get(mon, zero) - coeff
        if coeff:
            h[mon] = coeff
        elif mon in h:
            del h[mon]
    return h


def smp_degree(d: smp[Er], i_gen: int, n: int, domain: Domain[Er]) -> int:
    # Return the leading degree in one variable.
    if not d:
        return -1
    elif i_gen < 0:
        return 0
    else:
        return max(mon[i_gen] for mon in d)


def smp_degrees(d: smp[Er], n: int, domain: Domain[Er]) -> tuple[int, ...]:
    # Return the leading degrees in all variables.
    if not d:
        return (-1,) * n
    else:
        return tuple(max(mon[i] for mon in d) for i in range(n))


def smp_tail_degree(d: smp[Er], i_gen: int, n: int, domain: Domain[Er]) -> int:
    # Return the tail degree in one variable.
    if not d:
        return -1
    elif i_gen < 0:
        return 0
    else:
        return min(mon[i_gen] for mon in d)


def smp_tail_degrees(d: smp[Er], n: int, domain: Domain[Er]) -> tuple[int, ...]:
    # Return the tail degrees in all variables.
    if not d:
        return (-1,) * n
    else:
        return tuple(min(mon[i] for mon in d) for i in range(n))


def smp_diff(d: smp[Er], i_gen: int, n: int, domain: Domain[Er]) -> smp[Er]:
    # Differentiate a sparse polynomial with respect to one variable.
    h: smp[Er] = {}
    for mon, coeff in d.items():
        exp = mon[i_gen]
        if exp:
            new_mon = mon[:i_gen] + (exp - 1,) + mon[i_gen + 1 :]
            new_coeff = domain.convert(coeff * exp)
            if new_coeff:
                h[new_mon] = new_coeff
    return h


def smp_coeff_wrt(
    d: smp[Er], x_index: int, deg: int, n: int, domain: Domain[Er]
) -> smp[Er]:
    # Extract the coefficient wrt one variable of a fixed degree.
    h: smp[Er] = {}
    for mon, coeff in d.items():
        if mon[x_index] == deg:
            new_mon = mon[:x_index] + (0,) + mon[x_index + 1 :]
            h[new_mon] = coeff
    return h


def smp_deflate(
    polys: list[smp[Er]], n: int, domain: Domain[Er]
) -> tuple[tuple[int, ...], list[smp[Er]]]:
    # Deflate sparse polynomial exponents by their coordinate-wise gcd.
    J = [0] * n
    for d in polys:
        for mon in d:
            for i, exp in enumerate(mon):
                J[i] = igcd(J[i], exp)

    for i, gcd_exp in enumerate(J):
        if not gcd_exp:
            J[i] = 1

    J_tuple = tuple(J)

    if all(gcd_exp == 1 for gcd_exp in J_tuple):
        return J_tuple, polys

    return J_tuple, smp_deflate_by(polys, J_tuple, n, domain)


def smp_deflate_by(
    polys: list[smp[Er]], J: Sequence[int], n: int, domain: Domain[Er]
) -> list[smp[Er]]:
    # Deflate sparse polynomial exponents by given coordinate-wise factors.
    deflated = []
    for d in polys:
        h: smp[Er] = {}
        for mon, coeff in d.items():
            new_mon = tuple(exp // gcd_exp for exp, gcd_exp in zip(mon, J))
            h[new_mon] = coeff
        deflated.append(h)

    return deflated


def smp_inflate(d: smp[Er], J: Sequence[int], n: int, domain: Domain[Er]) -> smp[Er]:
    # Inflate sparse polynomial exponents by coordinate-wise factors.
    h: smp[Er] = {}
    for mon, coeff in d.items():
        new_mon = tuple(exp * factor for exp, factor in zip(mon, J))
        h[new_mon] = coeff
    return h


def smp_content(d: smp[Er], n: int, domain: Domain[Er]) -> Er:
    # Return the gcd of all sparse polynomial coefficients.
    cont = domain.zero
    gcd = domain.gcd

    for coeff in d.values():
        cont = gcd(cont, coeff)

    return cont


def smp_primitive(d: smp[Er], n: int, domain: Domain[Er]) -> tuple[Er, smp[Er]]:
    # Return the content and primitive part of a sparse polynomial.
    cont = smp_content(d, n, domain)
    if cont == domain.zero:
        return cont, d.copy()
    return cont, smp_quo_ground(d, cont, n, domain)


def smp_clear_denoms(d: smp[Er], n: int, domain: Domain[Er]) -> tuple[Er, smp[Er]]:
    # Clear denominators from sparse polynomial coefficients.
    if not domain.is_Field or not domain.has_assoc_Ring:
        return domain.one, d.copy()

    ground_ring = domain.get_ring()
    common = ground_ring.one
    lcm = ground_ring.lcm
    denom = domain.denom

    for coeff in d.values():
        common = lcm(common, denom(coeff))

    h: smp[Er] = {}
    for mon, coeff in d.items():
        new_coeff = coeff * common
        if new_coeff:
            h[mon] = new_coeff

    return common, h


def smp_is_zero(d: smp[Er], n: int, domain: Domain[Er]) -> bool:
    # Return whether a sparse polynomial is zero.
    return not d


def smp_is_one(d: smp[Er], n: int, domain: Domain[Er]) -> bool:
    # Return whether a sparse polynomial is one.
    zm = (0,) * n
    return len(d) == 1 and domain.is_one(d.get(zm, domain.zero))


def smp_is_ground(d: smp[Er], n: int, domain: Domain[Er]) -> bool:
    # Return whether a sparse polynomial is constant.
    return not d or (len(d) == 1 and (0,) * n in d)


def smp_is_term(d: smp[Er], n: int, domain: Domain[Er]) -> bool:
    # Return whether a sparse polynomial has at most one term.
    return len(d) <= 1


def smp_is_monomial(
    d: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex
) -> bool:
    # Return whether a sparse polynomial is zero or one monic term.
    return not d or (len(d) == 1 and smp_LC(d, n, domain, order) == 1)


def smp_is_negative(
    d: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex
) -> bool:
    # Return whether the leading coefficient is negative.
    return domain.is_negative(smp_LC(d, n, domain, order))


def smp_is_positive(
    d: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex
) -> bool:
    # Return whether the leading coefficient is positive.
    return domain.is_positive(smp_LC(d, n, domain, order))


def smp_is_nonnegative(
    d: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex
) -> bool:
    # Return whether the leading coefficient is nonnegative.
    return domain.is_nonnegative(smp_LC(d, n, domain, order))


def smp_is_nonpositive(
    d: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex
) -> bool:
    # Return whether the leading coefficient is nonpositive.
    return domain.is_nonpositive(smp_LC(d, n, domain, order))


def smp_is_monic(
    d: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex
) -> bool:
    # Return whether the leading coefficient is one.
    return domain.is_one(smp_LC(d, n, domain, order))


def smp_is_primitive(d: smp[Er], n: int, domain: Domain[Er]) -> bool:
    # Return whether the content is one.
    return domain.is_one(smp_content(d, n, domain))


def smp_is_linear(d: smp[Er], n: int, domain: Domain[Er]) -> bool:
    # Return whether every monomial has total degree at most one.
    return all(sum(mon) <= 1 for mon in d)


def smp_is_quadratic(d: smp[Er], n: int, domain: Domain[Er]) -> bool:
    # Return whether every monomial has total degree at most two.
    return all(sum(mon) <= 2 for mon in d)


def smp_iadd_monom(
    d: smp[Er], term: tuple[monom, Er], n: int, domain: Domain[Er]
) -> smp[Er]:
    # Add one term to a sparse polynomial in place.
    mon, coeff = term
    coeff = d.get(mon, domain.zero) + coeff
    if coeff:
        d[mon] = coeff
    elif mon in d:
        del d[mon]
    return d


def smp_iadd_poly_monom(
    d1: smp[Er], d2: smp[Er], term: tuple[monom, Er], n: int, domain: Domain[Er]
) -> smp[Er]:
    # Add d2 multiplied by one term to d1 in place.
    term_mon, term_coeff = term
    zero = domain.zero

    for mon, coeff in d2.items():
        p_mon = tuple(mon[i] + term_mon[i] for i in range(n))
        c = d1.get(p_mon, zero) + coeff * term_coeff
        if c:
            d1[p_mon] = c
        elif p_mon in d1:
            del d1[p_mon]

    return d1


def smp_add_ground(d: smp[Er], c: Er, n: int, domain: Domain[Er]) -> smp[Er]:
    # Add a ground coefficient to the constant term.
    h = d.copy()
    if not c:
        return h

    zm: monom = (0,) * n
    coeff = d.get(zm, domain.zero) + c
    if coeff:
        h[zm] = coeff
    elif zm in h:
        del h[zm]
    return h


def smp_sub_ground(d: smp[Er], c: Er, n: int, domain: Domain[Er]) -> smp[Er]:
    # Subtract a ground coefficient from the constant term.
    h = d.copy()
    if not c:
        return h

    zm: monom = (0,) * n
    coeff = d.get(zm, domain.zero) - c
    if coeff:
        h[zm] = coeff
    elif zm in h:
        del h[zm]
    return h


def smp_rsub_ground(d: smp[Er], c: Er, n: int, domain: Domain[Er]) -> smp[Er]:
    # Subtract the sparse polynomial from a ground coefficient.
    h = {mon: -coeff for mon, coeff in d.items() if coeff}
    return smp_add_ground(h, c, n, domain)


def smp_mul_ground(d: smp[Er], x: Er, n: int, domain: Domain[Er]) -> smp[Er]:
    # Multiply every coefficient by a ground element.
    if not x:
        return {}

    h: smp[Er] = {}
    for mon, coeff in d.items():
        coeff *= x
        if coeff:
            h[mon] = coeff
    return h


def smp_imul_num(d: smp[Er], c: Er, n: int, domain: Domain[Er]) -> smp[Er]:
    # Multiply every coefficient by a ground element in place.
    if not c:
        d.clear()
        return d

    for mon in list(d):
        coeff = d[mon] * c
        if coeff:
            d[mon] = coeff
        else:
            del d[mon]

    return d


def smp_quo_ground(d: smp[Er], x: Er, n: int, domain: Domain[Er]) -> smp[Er]:
    # Divide exactly divisible coefficients by a ground element.
    h: smp[Er] = {}

    if domain.is_Field:
        quo = domain.quo
        for mon, coeff in d.items():
            quotient = quo(coeff, x)
            if quotient:
                h[mon] = quotient
    else:
        # This is not valid for all domains (e.g. GF(p)).
        for mon, coeff in d.items():
            if not (coeff % x):  # type: ignore
                quotient = coeff // x  # type: ignore
                if quotient:
                    h[mon] = quotient

    return h


def smp_trunc_ground(d: smp[Er], p: Er, n: int, domain: Domain[Er]) -> smp[Er]:
    # Reduce coefficients modulo a ground element.
    h: smp[Er] = {}

    if domain.is_ZZ:
        p = domain.convert(p)
        for mon, coeff in d.items():
            coeff = domain.convert(coeff % p)  # type: ignore
            if coeff > p // 2:  # type: ignore
                coeff = domain.convert(coeff - p)
            if coeff:
                h[mon] = coeff
    else:
        for mon, coeff in d.items():
            coeff = domain.convert(coeff % p)  # type: ignore
            if coeff:
                h[mon] = coeff

    return h


def smp_subs(
    d: smp[Er], subs_dict: Mapping[int, Er], n: int, domain: Domain[Er]
) -> smp[Er]:
    # Substitute selected variables while preserving monomial dimensions.
    h: smp[Er] = {}

    for mon, coeff in d.items():
        new_coeff = coeff
        new_mon_list = list(mon)

        for i, value in subs_dict.items():
            exp = mon[i]
            if exp > 0:
                new_coeff *= value**exp
            new_mon_list[i] = 0

        if new_coeff:
            new_mon = tuple(new_mon_list)
            coeff = h.get(new_mon, domain.zero) + new_coeff
            if coeff:
                h[new_mon] = coeff
            elif new_mon in h:
                del h[new_mon]

    return h


def smp_subs_drop(
    d: smp[Er], subs_dict: Mapping[int, Er], n: int, domain: Domain[Er]
) -> smp[Er]:
    # Substitute selected variables and remove their monomial coordinates.
    h: smp[Er] = {}
    kept_indices = [i for i in range(n) if i not in subs_dict]

    for mon, coeff in d.items():
        new_coeff = coeff

        for i, value in subs_dict.items():
            exp = mon[i]
            if exp > 0:
                new_coeff *= value**exp

        if new_coeff:
            new_mon = tuple(mon[i] for i in kept_indices)
            coeff = h.get(new_mon, domain.zero) + new_coeff
            if coeff:
                h[new_mon] = coeff
            elif new_mon in h:
                del h[new_mon]

    return h


def smp_evaluate(
    d: smp[Er], eval_dict: Mapping[int, Er], n: int, domain: Domain[Er]
) -> Er:
    # Evaluate all variables and return an element of the ground domain.
    result = domain.zero

    for mon, coeff in d.items():
        mon_value = domain.one
        for i in range(n):
            exp = mon[i]
            if exp > 0:
                mon_value *= eval_dict[i] ** exp

        result += coeff * mon_value

    return result


def smp_leading_expv(
    d: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex
) -> monom | None:
    # Return the greatest exponent tuple for the selected monomial order.
    if not d:
        return None
    elif order is lex:
        return max(d)
    else:
        return max(d, key=order)


def smp_get_coeff(d: smp[Er], expv: monom | None, n: int, domain: Domain[Er]) -> Er:
    # Return the coefficient of an exponent tuple or zero when absent.
    if expv is None:
        return domain.zero
    return d.get(expv, domain.zero)


def smp_LC(d: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex) -> Er:
    # Return the leading coefficient for the selected monomial order.
    return smp_get_coeff(d, smp_leading_expv(d, n, domain, order), n, domain)


def smp_LM(d: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex) -> monom:
    # Return the leading exponent tuple, using the zero monomial for zero.
    expv = smp_leading_expv(d, n, domain, order)
    if expv is None:
        return (0,) * n
    return expv


def smp_LT(
    d: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex
) -> tuple[monom, Er]:
    # Return the leading exponent tuple and coefficient.
    expv = smp_leading_expv(d, n, domain, order)
    if expv is None:
        return ((0,) * n, domain.zero)
    return (expv, smp_get_coeff(d, expv, n, domain))


def smp_leading_monom(
    d: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex
) -> smp[Er]:
    # Return the leading monomial as a one-term sparse polynomial.
    expv = smp_leading_expv(d, n, domain, order)
    if expv:
        return {expv: domain.one}
    return {}


def smp_leading_term(
    d: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex
) -> smp[Er]:
    # Return the leading term as a one-term sparse polynomial.
    expv = smp_leading_expv(d, n, domain, order)
    if expv is None:
        return {}
    return {expv: d[expv]}


def smp_const(d: smp[Er], n: int, domain: Domain[Er]) -> Er:
    # Return the constant coefficient of a sparse polynomial.
    return smp_get_coeff(d, (0,) * n, n, domain)


def smp_coeff(d: smp[Er], monomial: smp[Er], n: int, domain: Domain[Er]) -> Er:
    # Return the coefficient of a monomial represented by a sparse dictionary.
    terms = list(monomial.items())
    if len(terms) == 1:
        mon, coeff = terms[0]
        if coeff == domain.one:
            return smp_get_coeff(d, mon, n, domain)
    raise ValueError(f"expected a monomial, got {monomial}")


def smp_term_div(
    term1: tuple[monom, Er], term2: tuple[monom, Er], n: int, domain: Domain[Er]
) -> tuple[monom, Er] | None:
    # Divide one polynomial term by another when the division is exact.
    mon1, coeff1 = term1
    mon2, coeff2 = term2
    zm: monom = (0,) * n

    if mon2 == zm:
        mon = mon1
    else:
        mon = tuple(mon1[i] - mon2[i] for i in range(n))
        if any(exp < 0 for exp in mon):
            return None

    if domain.is_Field:
        return mon, domain.quo(coeff1, coeff2)
    elif not (coeff1 % coeff2):  # type: ignore
        return mon, domain.quo(coeff1, coeff2)
    else:
        return None


def smp_quo_term(
    d: smp[Er], term: tuple[monom, Er], n: int, domain: Domain[Er]
) -> smp[Er]:
    # Divide every exactly divisible term by a fixed polynomial term.
    h: smp[Er] = {}
    for f_term in d.items():
        q_term = smp_term_div(f_term, term, n, domain)
        if q_term is not None:
            mon, coeff = q_term
            if coeff:
                h[mon] = coeff
    return h


def smp_div_list(
    d: smp[Er],
    gs: list[smp[Er]],
    n: int,
    domain: Domain[Er],
    order: MonomialOrder = lex,
) -> tuple[list[smp[Er]], smp[Er]]:
    # Divide a sparse polynomial by a list of sparse polynomials.
    if not all(gs):
        raise ZeroDivisionError("polynomial division")

    qs: list[smp[Er]] = [{} for _ in gs]
    lts = [smp_LT(g, n, domain, order) for g in gs]
    f = d.copy()
    r: smp[Er] = {}

    while f:
        divoccurred = False
        lm_f = smp_leading_expv(f, n, domain, order)
        if lm_f is None:
            break

        for i, (g, lt_g) in enumerate(zip(gs, lts)):
            q_term = smp_term_div((lm_f, f[lm_f]), lt_g, n, domain)
            if q_term is not None:
                q_mon, c = q_term
                smp_iadd_monom(qs[i], q_term, n, domain)
                smp_iadd_poly_monom(f, g, (q_mon, -c), n, domain)
                divoccurred = True
                break

        if not divoccurred:
            smp_iadd_monom(r, (lm_f, f[lm_f]), n, domain)
            del f[lm_f]

    return qs, r


def smp_rem_list(
    d: smp[Er],
    gs: list[smp[Er]],
    n: int,
    domain: Domain[Er],
    order: MonomialOrder = lex,
) -> smp[Er]:
    # Compute the remainder after division by a list of polynomials.
    if not all(gs):
        raise ZeroDivisionError("polynomial division")

    r: smp[Er] = {}
    lts = [smp_LT(g, n, domain, order) for g in gs]
    f = d.copy()

    while f:
        lm_f = smp_leading_expv(f, n, domain, order)
        if lm_f is None:
            break
        lt_f = (lm_f, f[lm_f])

        for g, lt_g in zip(gs, lts):
            q_term = smp_term_div(lt_f, lt_g, n, domain)
            if q_term is not None:
                q_mon, c = q_term
                for mon, coeff in g.items():
                    p_mon = tuple(mon[i] + q_mon[i] for i in range(n))
                    c1 = f.get(p_mon, domain.zero) - c * coeff
                    if c1:
                        f[p_mon] = c1
                    elif p_mon in f:
                        del f[p_mon]
                break
        else:
            smp_iadd_monom(r, lt_f, n, domain)
            del f[lm_f]

    return r


def smp_compose(
    d: smp[Er],
    reps: list[tuple[int, smp[Er]]],
    initial_poly: smp[Er],
    n: int,
    domain: Domain[Er],
) -> smp[Er]:
    # Compose a sparse polynomial with replacement polynomials.
    poly = initial_poly.copy()

    for mon, coeff in d.items():
        mon_list = list(mon)
        subpoly: smp[Er] = {(0,) * n: domain.one}

        for i, g in reps:
            exp, mon_list[i] = mon_list[i], 0
            if exp:
                subpoly = smp_mul(
                    subpoly, smp_pow_generic(g, exp, domain, n), domain, n
                )

        smp_iadd_poly_monom(poly, subpoly, (tuple(mon_list), coeff), n, domain)

    return poly


def _smp_symmetric_poly(k: int, n: int, domain: Domain[Er]) -> smp[Er]:
    # Build the elementary symmetric polynomial of degree k.
    if k < 0 or k > n:
        raise ValueError(f"Cannot generate symmetric polynomial of order {k}")
    elif not k:
        return {(0,) * n: domain.one}

    h: smp[Er] = {}
    for subset in combinations(range(n), k):
        mon = tuple(int(i in subset) for i in range(n))
        h[mon] = domain.one
    return h


def smp_symmetrize(
    d: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex
) -> tuple[smp[Er], smp[Er], list[smp[Er]]]:
    # Rewrite a symmetric polynomial in elementary symmetric polynomials.
    f = d.copy()

    if not n:
        return f, {}, []

    polys = [_smp_symmetric_poly(i + 1, n, domain) for i in range(n)]
    poly_powers: dict[tuple[int, int], smp[Er]] = {}

    def get_poly_power(i: int, exp: int) -> smp[Er]:
        if (i, exp) not in poly_powers:
            poly_powers[(i, exp)] = smp_pow_generic(polys[i], exp, domain, n)
        return poly_powers[(i, exp)]

    indices = list(range(n - 1))
    weights = list(range(n, 0, -1))

    symmetric: smp[Er] = {}

    while f:
        height = -1
        mc: tuple[monom, Er] | None = None

        if order is lex:
            terms = sorted(f.items(), key=lambda term: term[0], reverse=True)
        else:
            terms = sorted(f.items(), key=lambda term: order(term[0]), reverse=True)

        for mon, coeff in terms:
            if all(mon[i] >= mon[i + 1] for i in indices):
                mon_height = max(
                    weight * mon_exp for weight, mon_exp in zip(weights, mon)
                )

                if mon_height > height:
                    height = mon_height
                    mc = mon, coeff

        if mc is None:
            break
        mon, coeff = mc

        exps = []
        for m1, m2 in zip(mon, mon[1:] + (0,)):
            exps.append(m1 - m2)

        smp_iadd_monom(symmetric, (tuple(exps), coeff), n, domain)

        h: smp[Er] = {(0,) * n: coeff}
        for i, exp in enumerate(exps):
            h = smp_mul(h, get_poly_power(i, exp), domain, n)

        f = smp_sub(f, h, domain, n)

    return symmetric, f, polys


def smp_square(d: smp[Er], domain: Domain[Er], n: int) -> smp[Er]:
    # Square a sparse polynomial using cross and diagonal terms.
    zero = domain.zero
    h: smp[Er] = {}
    get = h.get
    mons = list(d.keys())
    for i, mon1 in enumerate(mons):
        coeff1 = d[mon1]
        for j in range(i):
            mon2 = mons[j]
            mon = tuple(mon1[i] + mon2[i] for i in range(n))
            coeff = get(mon, zero) + coeff1 * d[mon2]
            if coeff:
                h[mon] = coeff
            elif mon in h:
                del h[mon]

    for mon, coeff in list(h.items()):
        coeff *= 2
        if coeff:
            h[mon] = coeff
        else:
            del h[mon]

    get = h.get
    for mon, coeff in d.items():
        mon2 = tuple(2 * mon[i] for i in range(n))
        coeff = get(mon2, zero) + coeff**2
        if coeff:
            h[mon2] = coeff
        elif mon2 in h:
            del h[mon2]

    return h


def smp_pow_generic(d: smp[Er], exp: int, domain: Domain[Er], n: int) -> smp[Er]:
    # Raise a sparse polynomial to a power by binary exponentiation.
    if exp == 0:
        if not d:
            raise ValueError("0**0")
        um: monom = (0,) * n
        return {um: domain.one}
    elif exp < 0:
        raise ValueError(f"exponent must be a non-negative integer, got {exp}")
    if not d:
        return {}

    zm: monom = (0,) * n
    h: smp[Er] = {zm: domain.one}
    f = d

    while True:
        if exp & 1:
            h = smp_mul(h, f, domain, n)
            exp -= 1
            if not exp:
                break

        f = smp_square(f, domain, n)
        exp = exp // 2

    return h


def smp_pow_multinomial(d: smp[Er], exp: int, domain: Domain[Er], n: int) -> smp[Er]:
    # Raise a sparse polynomial to a power by multinomial expansion.
    if exp == 0:
        if not d:
            raise ValueError("0**0")
        um: monom = (0,) * n
        return {um: domain.one}
    elif exp < 0:
        raise ValueError(f"exponent must be a non-negative integer, got {exp}")

    if not d:
        return {}

    multinomials = multinomial_coefficients(len(d), exp).items()
    zm: monom = (0,) * n
    terms = d.items()
    zero = domain.zero
    h: smp[Er] = {}

    for multinomial, multinomial_coeff in multinomials:
        p_mon = zm
        c = domain.convert(multinomial_coeff)
        for term_exp, (mon, coeff) in zip(multinomial, terms):
            if term_exp:
                p_mon = tuple(p_mon[i] + mon[i] * term_exp for i in range(n))
                c *= coeff**term_exp

        coeff = h.get(p_mon, zero) + c

        if coeff:
            h[p_mon] = coeff
        elif p_mon in h:
            del h[p_mon]

    return h
