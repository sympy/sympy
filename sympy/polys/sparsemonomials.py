"""Sparse monomial helpers for the experimental SMP representation.

A sparse monomial is represented as a tuple of ``(generator_index, exponent)``
pairs with strictly increasing generator indices and positive exponents::

    x0**3*x8*x100**5  ->  ((0, 3), (8, 1), (100, 5))

The empty tuple ``()`` is the constant monomial.

This module is intentionally independent of ``sparsetools`` so that the
existing dense-exponent sparse-polynomial algorithms can remain unchanged.
"""

from __future__ import annotations

from typing import Mapping, TypeVar

from sympy.polys.monomials import monom

_T = TypeVar("_T")

smonom = tuple[tuple[int, int], ...]
ssmp = dict[smonom, _T]


def from_dense(mon: monom) -> smonom:
    """Convert a dense exponent tuple to a sparse monomial."""
    return tuple((i, exp) for i, exp in enumerate(mon) if exp)


def to_dense(mon: smonom, n: int) -> monom:
    """Convert a sparse monomial to a dense exponent tuple."""
    result = [0] * n
    for i, exp in mon:
        result[i] = exp
    return tuple(result)


def degree(mon: smonom, i_gen: int) -> int:
    """Return the exponent of one generator."""
    for i, exp in mon:
        if i == i_gen:
            return exp
        if i > i_gen:
            break
    return 0


def set_exp(mon: smonom, i_gen: int, exp: int) -> smonom:
    """Set one exponent, inserting or removing its sparse entry."""
    result = []
    inserted = False

    for i, old_exp in mon:
        if i < i_gen:
            result.append((i, old_exp))
        elif i == i_gen:
            if exp:
                result.append((i, exp))
            inserted = True
        else:
            if not inserted:
                if exp:
                    result.append((i_gen, exp))
                inserted = True
            result.append((i, old_exp))

    if not inserted and exp:
        result.append((i_gen, exp))

    return tuple(result)


def mul_monom(a: smonom, b: smonom) -> smonom:
    """Multiply two sparse monomials."""
    i = j = 0
    result = []

    while i < len(a) and j < len(b):
        ai, ae = a[i]
        bi, be = b[j]

        if ai < bi:
            result.append((ai, ae))
            i += 1
        elif bi < ai:
            result.append((bi, be))
            j += 1
        else:
            result.append((ai, ae + be))
            i += 1
            j += 1

    result.extend(a[i:])
    result.extend(b[j:])
    return tuple(result)


def _lex_gt(a: smonom, b: smonom) -> bool:
    """Compare sparse monomials in the dense lexicographic order."""
    i = j = 0

    while i < len(a) and j < len(b):
        ai, ae = a[i]
        bi, be = b[j]

        if ai == bi:
            if ae != be:
                return ae > be
            i += 1
            j += 1
        elif ai < bi:
            # a has a positive exponent where b has zero.
            return True
        else:
            # b has a positive exponent where a has zero.
            return False

    return i < len(a)


def add(f: ssmp[_T], g: ssmp[_T], domain) -> ssmp[_T]:
    h = f.copy()
    zero = domain.zero

    for mon, coeff in g.items():
        coeff = h.get(mon, zero) + coeff
        if coeff:
            h[mon] = coeff
        elif mon in h:
            del h[mon]

    return h


def sub(f: ssmp[_T], g: ssmp[_T], domain) -> ssmp[_T]:
    h = f.copy()
    zero = domain.zero

    for mon, coeff in g.items():
        coeff = h.get(mon, zero) - coeff
        if coeff:
            h[mon] = coeff
        elif mon in h:
            del h[mon]

    return h


def neg(f: ssmp[_T]) -> ssmp[_T]:
    return {mon: -coeff for mon, coeff in f.items() if coeff}


def add_ground(d: ssmp[_T], c, domain) -> ssmp[_T]:
    h = d.copy()
    coeff = h.get((), domain.zero) + c

    if coeff:
        h[()] = coeff
    elif () in h:
        del h[()]

    return h


def sub_ground(d: ssmp[_T], c, domain) -> ssmp[_T]:
    return add_ground(d, -c, domain)


def mul_ground(d: ssmp[_T], c) -> ssmp[_T]:
    if not c:
        return {}

    h = {}
    for mon, coeff in d.items():
        coeff = coeff * c
        if coeff:
            h[mon] = coeff
    return h


def mul(f: ssmp[_T], g: ssmp[_T], domain) -> ssmp[_T]:
    zero = domain.zero
    h = {}

    for mon1, coeff1 in f.items():
        for mon2, coeff2 in g.items():
            mon = mul_monom(mon1, mon2)
            coeff = h.get(mon, zero) + coeff1 * coeff2

            if coeff:
                h[mon] = coeff
            elif mon in h:
                del h[mon]

    return h


def square(d: ssmp[_T], domain) -> ssmp[_T]:
    zero = domain.zero
    h = {}
    mons = list(d)

    # Cross terms.
    for i, mon1 in enumerate(mons):
        coeff1 = d[mon1]

        for j in range(i):
            mon = mul_monom(mon1, mons[j])
            coeff = h.get(mon, zero) + coeff1 * d[mons[j]]

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

    # Diagonal terms.
    for mon, coeff in d.items():
        mon2 = mul_monom(mon, mon)
        coeff2 = h.get(mon2, zero) + coeff**2

        if coeff2:
            h[mon2] = coeff2
        elif mon2 in h:
            del h[mon2]

    return h


def pow_generic(d: ssmp[_T], exp: int, domain) -> ssmp[_T]:
    if exp == 0:
        if not d:
            raise ValueError("0**0")
        return {(): domain.one}
    elif exp < 0:
        raise ValueError(f"exponent must be a non-negative integer, got {exp}")

    if not d:
        return {}

    h = {(): domain.one}
    f = d

    while True:
        if exp & 1:
            h = mul(h, f, domain)
            exp -= 1
            if not exp:
                break

        f = square(f, domain)
        exp //= 2

    return h


def poly_degree(d: ssmp[_T], i_gen: int) -> int:
    if not d:
        return -1
    elif i_gen < 0:
        return 0

    return max(degree(mon, i_gen) for mon in d)


def poly_degrees(d: ssmp[_T], n: int) -> tuple[int, ...]:
    if not d:
        return (-1,) * n

    result = [0] * n

    for mon in d:
        for i, exp in mon:
            if exp > result[i]:
                result[i] = exp

    return tuple(result)


def total_degree(d: ssmp[_T]) -> int:
    if not d:
        return -1

    return max(sum(exp for _, exp in mon) for mon in d)


def diff(d: ssmp[_T], i_gen: int, domain) -> ssmp[_T]:
    h = {}

    for mon, coeff in d.items():
        exp = degree(mon, i_gen)

        if exp:
            new_mon = set_exp(mon, i_gen, exp - 1)
            new_coeff = domain.convert(coeff * exp)

            if new_coeff:
                h[new_mon] = new_coeff

    return h


def integrate(d: ssmp[_T], i_gen: int, domain) -> ssmp[_T]:
    h = {}

    for mon, coeff in d.items():
        exp = degree(mon, i_gen) + 1
        new_mon = set_exp(mon, i_gen, exp)
        new_coeff = domain.quo(coeff, domain.convert(exp))

        if new_coeff:
            h[new_mon] = new_coeff

    return h


def LC(d: ssmp[_T], domain):
    """Return the leading coefficient in lex order."""
    if not d:
        return domain.zero

    lead = next(iter(d))

    for mon in d:
        if _lex_gt(mon, lead):
            lead = mon

    return d[lead]


def content(d: ssmp[_T], domain):
    cont = domain.zero

    for coeff in d.values():
        cont = domain.gcd(cont, coeff)

    return cont


def primitive(d: ssmp[_T], domain):
    cont = content(d, domain)

    if cont == domain.zero:
        return cont, d.copy()

    h = {}

    for mon, coeff in d.items():
        coeff = domain.exquo(coeff, cont)
        if coeff:
            h[mon] = coeff

    return cont, h


def clear_denoms(d: ssmp[_T], domain):
    if not domain.is_Field or not domain.has_assoc_Ring:
        return domain.one, d.copy()

    ground_ring = domain.get_ring()
    common = ground_ring.one

    for coeff in d.values():
        common = ground_ring.lcm(common, domain.denom(coeff))

    h = {}

    for mon, coeff in d.items():
        coeff = coeff * common
        if coeff:
            h[mon] = coeff

    return common, h


def trunc_ground(d: ssmp[_T], p, domain) -> ssmp[_T]:
    h = {}

    if domain.is_ZZ:
        t = domain.quo(p, domain.convert(2))

        for mon, coeff in d.items():
            coeff = domain.rem(coeff, p)

            if domain.is_positive(coeff - t):
                coeff -= p

            if coeff:
                h[mon] = coeff
    else:
        for mon, coeff in d.items():
            coeff = domain.rem(coeff, p)

            if coeff:
                h[mon] = coeff

    return h


def subs_drop(
    d: ssmp[_T], subs_dict: Mapping[int, _T], n: int, domain
) -> ssmp[_T]:
    """Substitute selected generators and remove those coordinates."""
    h = {}
    dropped = sorted(subs_dict)

    for mon, coeff in d.items():
        new_coeff = coeff
        new_mon = []

        for i, exp in mon:
            if i in subs_dict:
                new_coeff *= subs_dict[i] ** exp
            else:
                shift = 0
                for j in dropped:
                    if j < i:
                        shift += 1
                    else:
                        break
                new_mon.append((i - shift, exp))

        if new_coeff:
            new_mon = tuple(new_mon)
            coeff2 = h.get(new_mon, domain.zero) + new_coeff

            if coeff2:
                h[new_mon] = coeff2
            elif new_mon in h:
                del h[new_mon]

    return h


def is_zero(d: ssmp[_T]) -> bool:
    return not d


def is_one(d: ssmp[_T], domain) -> bool:
    return len(d) == 1 and domain.is_one(d.get((), domain.zero))


def is_ground(d: ssmp[_T]) -> bool:
    return not d or (len(d) == 1 and () in d)


def is_monic(d: ssmp[_T], domain) -> bool:
    return domain.is_one(LC(d, domain))


def is_primitive(d: ssmp[_T], domain) -> bool:
    return domain.is_one(content(d, domain))


def is_linear(d: ssmp[_T]) -> bool:
    return all(sum(exp for _, exp in mon) <= 1 for mon in d)


def is_quadratic(d: ssmp[_T]) -> bool:
    return all(sum(exp for _, exp in mon) <= 2 for mon in d)
