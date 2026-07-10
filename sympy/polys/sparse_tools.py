"""Low-level sparse polynomial operations on raw dict representations.

These helpers deliberately do not depend on PolyRing or PolyElement.
They operate on dicts mapping exponent tuples to coefficients and use
only the provided domain and number of variables.
"""

from __future__ import annotations

from itertools import combinations
from typing import TYPE_CHECKING, Mapping, Sequence

from sympy.core.intfunc import igcd
from sympy.ntheory.multinomial import multinomial_coefficients
from sympy.polys.orderings import lex, MonomialOrder

if TYPE_CHECKING:
    from sympy.polys.domains.domain import Domain, Er

ninf = float("-inf")


def dict_add(f: dict[tuple[int, ...], Er], g: dict[tuple[int, ...], Er],
    domain: Domain[Er], n: int) -> dict[tuple[int, ...], Er]:
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


def dict_mul(f: dict[tuple[int, ...], Er], g: dict[tuple[int, ...], Er],
    domain: Domain[Er], n: int) -> dict[tuple[int, ...], Er]:
    # Multiply two sparse polynomials and return a new dictionary.
    zero = domain.zero
    h: dict[tuple[int, ...], Er] = {}
    for mon1, coeff1 in f.items():
        for mon2, coeff2 in g.items():
            mon = tuple(mon1[i] + mon2[i] for i in range(n))
            coeff = h.get(mon, zero) + coeff1*coeff2
            if coeff:
                h[mon] = coeff
            elif mon in h:
                del h[mon]
    return h


def dict_sub(f: dict[tuple[int, ...], Er], g: dict[tuple[int, ...], Er],
    domain: Domain[Er], n: int) -> dict[tuple[int, ...], Er]:
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


def dict_degree(d: dict[tuple[int, ...], Er], i_gen: int, n: int,
    domain: Domain[Er]) -> float:
    # Return the leading degree in one variable.
    if not d:
        return ninf
    elif i_gen < 0:
        return 0
    else:
        return max(monom[i_gen] for monom in d)


def dict_degrees(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> tuple[float, ...]:
    # Return the leading degrees in all variables.
    if not d:
        return (ninf,)*n
    else:
        return tuple(max(monom[i] for monom in d) for i in range(n))


def dict_tail_degree(d: dict[tuple[int, ...], Er], i_gen: int, n: int,
    domain: Domain[Er]) -> float:
    # Return the tail degree in one variable.
    if not d:
        return ninf
    elif i_gen < 0:
        return 0
    else:
        return min(monom[i_gen] for monom in d)


def dict_tail_degrees(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> tuple[float, ...]:
    # Return the tail degrees in all variables.
    if not d:
        return (ninf,)*n
    else:
        return tuple(min(monom[i] for monom in d) for i in range(n))


def dict_diff(d: dict[tuple[int, ...], Er], i_gen: int, n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Differentiate a sparse polynomial with respect to one variable.
    h: dict[tuple[int, ...], Er] = {}
    for mon, coeff in d.items():
        exp = mon[i_gen]
        if exp:
            new_mon = mon[:i_gen] + (exp - 1,) + mon[i_gen + 1:]
            new_coeff = domain.convert(coeff*exp)
            if new_coeff:
                h[new_mon] = new_coeff
    return h


def dict_coeff_wrt(d: dict[tuple[int, ...], Er], x_index: int,
    deg: int, n: int, domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Extract the coefficient wrt one variable of a fixed degree.
    h: dict[tuple[int, ...], Er] = {}
    for mon, coeff in d.items():
        if mon[x_index] == deg:
            new_mon = mon[:x_index] + (0,) + mon[x_index + 1:]
            h[new_mon] = coeff
    return h


def dict_deflate(dicts_list: list[dict[tuple[int, ...], Er]], n: int,
    domain: Domain[Er]) -> tuple[tuple[int, ...],
        list[dict[tuple[int, ...], Er]]]:
    # Deflate sparse polynomial exponents by their coordinate-wise gcd.
    J = [0]*n
    for d in dicts_list:
        for mon in d:
            for i, exp in enumerate(mon):
                J[i] = igcd(J[i], exp)

    for i, gcd_exp in enumerate(J):
        if not gcd_exp:
            J[i] = 1

    J_tuple = tuple(J)

    if all(gcd_exp == 1 for gcd_exp in J_tuple):
        return J_tuple, dicts_list

    return J_tuple, dict_deflate_by(dicts_list, J_tuple, n, domain)


def dict_deflate_by(dicts_list: list[dict[tuple[int, ...], Er]],
    J: Sequence[int], n: int,
    domain: Domain[Er]) -> list[dict[tuple[int, ...], Er]]:
    # Deflate sparse polynomial exponents by given coordinate-wise factors.
    deflated = []
    for d in dicts_list:
        h: dict[tuple[int, ...], Er] = {}
        for mon, coeff in d.items():
            new_mon = tuple(exp // gcd_exp
                for exp, gcd_exp in zip(mon, J))
            h[new_mon] = coeff
        deflated.append(h)

    return deflated


def dict_inflate(d: dict[tuple[int, ...], Er], J: Sequence[int], n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Inflate sparse polynomial exponents by coordinate-wise factors.
    h: dict[tuple[int, ...], Er] = {}
    for mon, coeff in d.items():
        new_mon = tuple(exp*factor for exp, factor in zip(mon, J))
        h[new_mon] = coeff
    return h


def dict_content(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> Er:
    # Return the gcd of all sparse polynomial coefficients.
    cont = domain.zero
    gcd = domain.gcd

    for coeff in d.values():
        cont = gcd(cont, coeff)

    return cont


def dict_primitive(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> tuple[Er, dict[tuple[int, ...], Er]]:
    # Return the content and primitive part of a sparse polynomial.
    cont = dict_content(d, n, domain)
    if cont == domain.zero:
        return cont, d.copy()
    return cont, dict_quo_ground(d, cont, n, domain)


def dict_clear_denoms(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> tuple[Er, dict[tuple[int, ...], Er]]:
    # Clear denominators from sparse polynomial coefficients.
    if not domain.is_Field or not domain.has_assoc_Ring:
        return domain.one, d.copy()

    ground_ring = domain.get_ring()
    common = ground_ring.one
    lcm = ground_ring.lcm
    denom = domain.denom

    for coeff in d.values():
        common = lcm(common, denom(coeff))

    h: dict[tuple[int, ...], Er] = {}
    for mon, coeff in d.items():
        new_coeff = coeff*common
        if new_coeff:
            h[mon] = new_coeff

    return common, h


def dict_is_zero(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> bool:
    # Return whether a sparse polynomial is zero.
    return not d


def dict_is_one(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> bool:
    # Return whether a sparse polynomial is one.
    zero_monom = (0,)*n
    return len(d) == 1 and domain.is_one(d.get(zero_monom, domain.zero))


def dict_is_ground(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> bool:
    # Return whether a sparse polynomial is constant.
    return not d or (len(d) == 1 and (0,)*n in d)


def dict_is_term(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> bool:
    # Return whether a sparse polynomial has at most one term.
    return len(d) <= 1


def dict_is_monomial(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er], order: MonomialOrder = lex) -> bool:
    # Return whether a sparse polynomial is zero or one monic term.
    return not d or (len(d) == 1 and dict_LC(d, n, domain, order) == 1)


def dict_is_negative(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er], order: MonomialOrder = lex) -> bool:
    # Return whether the leading coefficient is negative.
    return domain.is_negative(dict_LC(d, n, domain, order))


def dict_is_positive(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er], order: MonomialOrder = lex) -> bool:
    # Return whether the leading coefficient is positive.
    return domain.is_positive(dict_LC(d, n, domain, order))


def dict_is_nonnegative(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er], order: MonomialOrder = lex) -> bool:
    # Return whether the leading coefficient is nonnegative.
    return domain.is_nonnegative(dict_LC(d, n, domain, order))


def dict_is_nonpositive(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er], order: MonomialOrder = lex) -> bool:
    # Return whether the leading coefficient is nonpositive.
    return domain.is_nonpositive(dict_LC(d, n, domain, order))


def dict_is_monic(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er], order: MonomialOrder = lex) -> bool:
    # Return whether the leading coefficient is one.
    return domain.is_one(dict_LC(d, n, domain, order))


def dict_is_primitive(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> bool:
    # Return whether the content is one.
    return domain.is_one(dict_content(d, n, domain))


def dict_is_linear(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> bool:
    # Return whether every monomial has total degree at most one.
    return all(sum(mon) <= 1 for mon in d)


def dict_is_quadratic(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> bool:
    # Return whether every monomial has total degree at most two.
    return all(sum(mon) <= 2 for mon in d)


def dict_iadd_monom(d: dict[tuple[int, ...], Er],
    term: tuple[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Add one term to a sparse polynomial in place.
    mon, coeff = term
    coeff = d.get(mon, domain.zero) + coeff
    if coeff:
        d[mon] = coeff
    elif mon in d:
        del d[mon]
    return d


def dict_iadd_poly_monom(d1: dict[tuple[int, ...], Er],
    d2: dict[tuple[int, ...], Er],
    term: tuple[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Add d2 multiplied by one term to d1 in place.
    term_mon, term_coeff = term
    zero = domain.zero

    for mon, coeff in d2.items():
        product_mon = tuple(mon[i] + term_mon[i] for i in range(n))
        product_coeff = d1.get(product_mon, zero) + coeff*term_coeff
        if product_coeff:
            d1[product_mon] = product_coeff
        elif product_mon in d1:
            del d1[product_mon]

    return d1


def dict_add_ground(d: dict[tuple[int, ...], Er], c: Er, n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Add a ground coefficient to the constant term.
    h = d.copy()
    if not c:
        return h

    zero_monom: tuple[int, ...] = (0,)*n
    coeff = d.get(zero_monom, domain.zero) + c
    if coeff:
        h[zero_monom] = coeff
    elif zero_monom in h:
        del h[zero_monom]
    return h



def dict_sub_ground(d: dict[tuple[int, ...], Er], c: Er, n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Subtract a ground coefficient from the constant term.
    h = d.copy()
    if not c:
        return h

    zero_monom: tuple[int, ...] = (0,)*n
    coeff = d.get(zero_monom, domain.zero) - c
    if coeff:
        h[zero_monom] = coeff
    elif zero_monom in h:
        del h[zero_monom]
    return h


def dict_rsub_ground(d: dict[tuple[int, ...], Er], c: Er, n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Subtract the sparse polynomial from a ground coefficient.
    h = {mon: -coeff for mon, coeff in d.items() if coeff}
    return dict_add_ground(h, c, n, domain)


def dict_mul_ground(d: dict[tuple[int, ...], Er], x: Er, n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Multiply every coefficient by a ground element.
    if not x:
        return {}

    h: dict[tuple[int, ...], Er] = {}
    for mon, coeff in d.items():
        product = coeff*x
        if product:
            h[mon] = product
    return h


def dict_imul_num(d: dict[tuple[int, ...], Er], c: Er, n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Multiply every coefficient by a ground element in place.
    if not c:
        d.clear()
        return d

    for mon in list(d):
        coeff = d[mon]*c
        if coeff:
            d[mon] = coeff
        else:
            del d[mon]

    return d


def dict_quo_ground(d: dict[tuple[int, ...], Er], x: Er, n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Divide exactly divisible coefficients by a ground element.
    h: dict[tuple[int, ...], Er] = {}

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


def dict_trunc_ground(d: dict[tuple[int, ...], Er], p: Er, n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Reduce coefficients modulo a ground element.
    h: dict[tuple[int, ...], Er] = {}

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


def dict_subs(d: dict[tuple[int, ...], Er],
    subs_dict: Mapping[int, Er], n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Substitute selected variables while preserving monomial dimensions.
    h: dict[tuple[int, ...], Er] = {}

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


def dict_subs_drop(d: dict[tuple[int, ...], Er],
    subs_dict: Mapping[int, Er], n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Substitute selected variables and remove their monomial coordinates.
    h: dict[tuple[int, ...], Er] = {}
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


def dict_evaluate(d: dict[tuple[int, ...], Er],
    eval_dict: Mapping[int, Er], n: int, domain: Domain[Er]) -> Er:
    # Evaluate all variables and return an element of the ground domain.
    result = domain.zero

    for mon, coeff in d.items():
        mon_value = domain.one
        for i in range(n):
            exp = mon[i]
            if exp > 0:
                mon_value *= eval_dict[i]**exp

        result += coeff*mon_value

    return result


def dict_leading_expv(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er],
    order: MonomialOrder = lex) -> tuple[int, ...] | None:
    # Return the greatest exponent tuple for the selected monomial order.
    if not d:
        return None
    elif order is lex:
        return max(d)
    else:
        return max(d, key=order)


def dict_get_coeff(d: dict[tuple[int, ...], Er],
    expv: tuple[int, ...] | None, n: int, domain: Domain[Er]) -> Er:
    # Return the coefficient of an exponent tuple or zero when absent.
    if expv is None:
        return domain.zero
    return d.get(expv, domain.zero)


def dict_LC(d: dict[tuple[int, ...], Er], n: int, domain: Domain[Er],
    order: MonomialOrder = lex) -> Er:
    # Return the leading coefficient for the selected monomial order.
    return dict_get_coeff(
        d, dict_leading_expv(d, n, domain, order), n, domain)


def dict_LM(d: dict[tuple[int, ...], Er], n: int, domain: Domain[Er],
    order: MonomialOrder = lex) -> tuple[int, ...]:
    # Return the leading exponent tuple, using the zero monomial for zero.
    expv = dict_leading_expv(d, n, domain, order)
    if expv is None:
        return (0,)*n
    return expv


def dict_LT(d: dict[tuple[int, ...], Er], n: int, domain: Domain[Er],
    order: MonomialOrder = lex) -> tuple[tuple[int, ...], Er]:
    # Return the leading exponent tuple and coefficient.
    expv = dict_leading_expv(d, n, domain, order)
    if expv is None:
        return ((0,)*n, domain.zero)
    return (expv, dict_get_coeff(d, expv, n, domain))


def dict_leading_monom(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er],
    order: MonomialOrder = lex) -> dict[tuple[int, ...], Er]:
    # Return the leading monomial as a one-term sparse polynomial.
    expv = dict_leading_expv(d, n, domain, order)
    if expv:
        return {expv: domain.one}
    return {}


def dict_leading_term(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er],
    order: MonomialOrder = lex) -> dict[tuple[int, ...], Er]:
    # Return the leading term as a one-term sparse polynomial.
    expv = dict_leading_expv(d, n, domain, order)
    if expv is None:
        return {}
    return {expv: d[expv]}


def dict_const(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> Er:
    # Return the constant coefficient of a sparse polynomial.
    return dict_get_coeff(d, (0,)*n, n, domain)


def dict_coeff(d: dict[tuple[int, ...], Er],
    element_dict: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> Er:
    # Return the coefficient of a monomial represented by a sparse dictionary.
    terms = list(element_dict.items())
    if len(terms) == 1:
        mon, coeff = terms[0]
        if coeff == domain.one:
            return dict_get_coeff(d, mon, n, domain)
    raise ValueError(f"expected a monomial, got {element_dict}")


def dict_term_div(term1: tuple[tuple[int, ...], Er],
    term2: tuple[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> tuple[tuple[int, ...], Er] | None:
    # Divide one polynomial term by another when the division is exact.
    mon1, coeff1 = term1
    mon2, coeff2 = term2
    zero_monom: tuple[int, ...] = (0,)*n

    if mon2 == zero_monom:
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


def dict_quo_term(d: dict[tuple[int, ...], Er],
    term: tuple[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Divide every exactly divisible term by a fixed polynomial term.
    h: dict[tuple[int, ...], Er] = {}
    for dividend_term in d.items():
        quotient_term = dict_term_div(
            dividend_term, term, n, domain)
        if quotient_term is not None:
            mon, coeff = quotient_term
            if coeff:
                h[mon] = coeff
    return h


def dict_div_list(d: dict[tuple[int, ...], Er],
    divisors: list[dict[tuple[int, ...], Er]], n: int,
    domain: Domain[Er],
    order: MonomialOrder = lex
    ) -> tuple[list[dict[tuple[int, ...], Er]],
               dict[tuple[int, ...], Er]]:
    # Divide a sparse polynomial by a list of sparse polynomials.
    if not all(divisors):
        raise ZeroDivisionError("polynomial division")

    quotients: list[dict[tuple[int, ...], Er]] = [
        {} for _ in divisors]
    divisor_terms = [
        dict_LT(divisor, n, domain, order)
        for divisor in divisors
    ]
    dividend = d.copy()
    remainder: dict[tuple[int, ...], Er] = {}

    while dividend:
        division_occurred = False
        dividend_expv = dict_leading_expv(
            dividend, n, domain, order)
        if dividend_expv is None:
            break

        for i, (divisor, divisor_term) in enumerate(
            zip(divisors, divisor_terms)):
            quotient_term = dict_term_div(
                (dividend_expv, dividend[dividend_expv]),
                divisor_term, n, domain)
            if quotient_term is not None:
                quotient_expv, quotient_coeff = quotient_term
                dict_iadd_monom(
                    quotients[i], quotient_term, n, domain)
                dict_iadd_poly_monom(
                    dividend, divisor,
                    (quotient_expv, -quotient_coeff), n, domain)
                division_occurred = True
                break

        if not division_occurred:
            dict_iadd_monom(
                remainder,
                (dividend_expv, dividend[dividend_expv]), n, domain)
            del dividend[dividend_expv]

    return quotients, remainder


def dict_rem_list(d: dict[tuple[int, ...], Er],
    divisors: list[dict[tuple[int, ...], Er]], n: int,
    domain: Domain[Er],
    order: MonomialOrder = lex) -> dict[tuple[int, ...], Er]:
    # Compute the remainder after division by a list of polynomials.
    if not all(divisors):
        raise ZeroDivisionError("polynomial division")

    remainder: dict[tuple[int, ...], Er] = {}
    divisor_terms = [
        dict_LT(divisor, n, domain, order)
        for divisor in divisors
    ]
    dividend = d.copy()

    while dividend:
        dividend_expv = dict_leading_expv(
            dividend, n, domain, order)
        if dividend_expv is None:
            break
        dividend_term = (dividend_expv, dividend[dividend_expv])

        for divisor, divisor_term in zip(divisors, divisor_terms):
            quotient_term = dict_term_div(
                dividend_term, divisor_term, n, domain)
            if quotient_term is not None:
                quotient_mon, quotient_coeff = quotient_term
                for mon, coeff in divisor.items():
                    product_mon = tuple(
                        mon[i] + quotient_mon[i] for i in range(n))
                    product_coeff = (
                        dividend.get(product_mon, domain.zero)
                        - quotient_coeff*coeff
                    )
                    if product_coeff:
                        dividend[product_mon] = product_coeff
                    elif product_mon in dividend:
                        del dividend[product_mon]
                break
        else:
            dict_iadd_monom(
                remainder, dividend_term, n, domain)
            del dividend[dividend_expv]

    return remainder


def dict_compose(d: dict[tuple[int, ...], Er],
    replacements: list[tuple[int, dict[tuple[int, ...], Er]]],
    initial_poly_dict: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Compose a sparse polynomial with replacement polynomials.
    poly = initial_poly_dict.copy()

    for monom, coeff in d.items():
        monom_list = list(monom)
        subpoly: dict[tuple[int, ...], Er] = {
            (0,)*n: domain.one
        }

        for i, g in replacements:
            exp, monom_list[i] = monom_list[i], 0
            if exp:
                subpoly = dict_mul(
                    subpoly, dict_pow_generic(g, exp, domain, n),
                    domain, n)

        dict_iadd_poly_monom(
            poly, subpoly, (tuple(monom_list), coeff), n, domain)

    return poly


def _dict_symmetric_poly(k: int, n: int,
    domain: Domain[Er]) -> dict[tuple[int, ...], Er]:
    # Build the elementary symmetric polynomial of degree k.
    if k < 0 or k > n:
        raise ValueError(
            f"Cannot generate symmetric polynomial of order {k}")
    elif not k:
        return {(0,)*n: domain.one}

    h: dict[tuple[int, ...], Er] = {}
    for subset in combinations(range(n), k):
        monom = tuple(int(i in subset) for i in range(n))
        h[monom] = domain.one
    return h


def dict_symmetrize(d: dict[tuple[int, ...], Er], n: int,
    domain: Domain[Er], order: MonomialOrder = lex) -> tuple[
        dict[tuple[int, ...], Er],
        dict[tuple[int, ...], Er],
        list[dict[tuple[int, ...], Er]]
    ]:
    # Rewrite a symmetric polynomial in elementary symmetric polynomials.
    f = d.copy()

    if not n:
        return f, {}, []

    polys = [_dict_symmetric_poly(i + 1, n, domain) for i in range(n)]
    poly_powers: dict[
        tuple[int, int], dict[tuple[int, ...], Er]
    ] = {}

    def get_poly_power(i: int, exp: int) -> dict[tuple[int, ...], Er]:
        if (i, exp) not in poly_powers:
            poly_powers[(i, exp)] = dict_pow_generic(
                polys[i], exp, domain, n)
        return poly_powers[(i, exp)]

    indices = list(range(n - 1))
    weights = list(range(n, 0, -1))

    symmetric: dict[tuple[int, ...], Er] = {}

    while f:
        height = -1
        leading_monom: tuple[int, ...] | None = None
        leading_coeff: Er | None = None

        if order is lex:
            terms = sorted(
                f.items(), key=lambda term: term[0], reverse=True)
        else:
            terms = sorted(
                f.items(), key=lambda term: order(term[0]), reverse=True)

        for monom, coeff in terms:
            if all(monom[i] >= monom[i + 1] for i in indices):
                monom_height = max(
                    weight*monom_exp
                    for weight, monom_exp in zip(weights, monom))

                if monom_height > height:
                    height = monom_height
                    leading_monom, leading_coeff = monom, coeff

        if height == -1:
            break

        assert leading_monom is not None
        assert leading_coeff is not None
        monom = leading_monom
        coeff = leading_coeff

        exponents = []
        for m1, m2 in zip(monom, monom[1:] + (0,)):
            exponents.append(m1 - m2)

        dict_iadd_monom(
            symmetric, (tuple(exponents), coeff), n, domain)

        product: dict[tuple[int, ...], Er] = {
            (0,)*n: coeff
        }
        for i, exp in enumerate(exponents):
            product = dict_mul(
                product, get_poly_power(i, exp), domain, n)

        f = dict_sub(f, product, domain, n)

    return symmetric, f, polys


def dict_square(d: dict[tuple[int, ...], Er], domain: Domain[Er],
    n: int) -> dict[tuple[int, ...], Er]:
    # Square a sparse polynomial using cross and diagonal terms.
    zero = domain.zero
    h: dict[tuple[int, ...], Er] = {}
    get = h.get
    mons = list(d.keys())
    for i, mon1 in enumerate(mons):
        coeff1 = d[mon1]
        for j in range(i):
            mon2 = mons[j]
            mon = tuple(mon1[i] + mon2[i] for i in range(n))
            coeff = get(mon, zero) + coeff1*d[mon2]
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
        mon2 = tuple(2*mon[i] for i in range(n))
        coeff = get(mon2, zero) + coeff**2
        if coeff:
            h[mon2] = coeff
        elif mon2 in h:
            del h[mon2]

    return h


def dict_pow_generic(d: dict[tuple[int, ...], Er], exp: int, domain: Domain[Er],
    n: int) -> dict[tuple[int, ...], Er]:
    # Raise a sparse polynomial to a power by binary exponentiation.
    if exp == 0:
        if not d:
            raise ValueError("0**0")
        unit_monom: tuple[int, ...] = (0,)*n
        return {unit_monom: domain.one}
    elif exp < 0:
        raise ValueError(
            f"exponent must be a non-negative integer, got {exp}")
    if not d:
        return {}

    zero_monom: tuple[int, ...] = (0,)*n
    p: dict[tuple[int, ...], Er] = {zero_monom: domain.one}
    c = d

    while True:
        if exp & 1:
            p = dict_mul(p, c, domain, n)
            exp -= 1
            if not exp:
                break

        c = dict_square(c, domain, n)
        exp = exp // 2

    return p


def dict_pow_multinomial(d: dict[tuple[int, ...], Er], exp: int,
    domain: Domain[Er], n: int) -> dict[tuple[int, ...], Er]:
    # Raise a sparse polynomial to a power by multinomial expansion.
    if exp == 0:
        if not d:
            raise ValueError("0**0")
        unit_monom: tuple[int, ...] = (0,)*n
        return {unit_monom: domain.one}
    elif exp < 0:
        raise ValueError(
            f"exponent must be a non-negative integer, got {exp}")

    if not d:
        return {}

    multinomials = multinomial_coefficients(len(d), exp).items()
    zero_monom: tuple[int, ...] = (0,)*n
    terms = d.items()
    zero = domain.zero
    h: dict[tuple[int, ...], Er] = {}

    for multinomial, multinomial_coeff in multinomials:
        product_monom = zero_monom
        product_coeff = domain.convert(multinomial_coeff)
        for term_exp, (mon, coeff) in zip(multinomial, terms):
            if term_exp:
                product_monom = tuple(
                    product_monom[i] + mon[i]*term_exp for i in range(n))
                product_coeff *= coeff**term_exp

        coeff = h.get(product_monom, zero) + product_coeff

        if coeff:
            h[product_monom] = coeff
        elif product_monom in h:
            del h[product_monom]

    return h
