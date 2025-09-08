from __future__ import annotations

from typing import Sequence, TYPE_CHECKING

if TYPE_CHECKING:
    from typing import TypeAlias, Union


from sympy.polys.densearith import (
    dup_add,
    dup_add_ground,
    dup_mul,
    dup_mul_ground,
    dup_neg,
    dup_sub,
    dup_sub_ground,
    dup_exquo,
    dup_series_mul,
    dup_series_sqr,
    dup_series_pow,
    dup_rshift,
)
from sympy.polys.densebasic import dup, dup_degree, dup_reverse, dup_truncate
from sympy.polys.densetools import (
    dup_diff,
    dup_integrate,
    dup_revert,
    dup_series_compose,
    dup_series_reversion,
)
from sympy.polys.polyerrors import NotReversible, ExactQuotientFailed
from sympy.polys.domains import Domain, QQ, ZZ
from sympy.polys.domains.domain import Er, Ef
from sympy.polys.domains.field import Field
from sympy.polys.series.base import series_pprint
from sympy.external.gmpy import MPZ, MPQ
from sympy.polys.ring_series import _giant_steps


USeries: TypeAlias = "tuple[list[Er], Union[int, None]]"


def _useries(
    coeffs: dup[Er], series_prec: int | None, dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    """Helper function to decide if the polynomial can be exact or it become a useries
    element."""

    deg = dup_degree(coeffs)
    if series_prec is None:
        if deg < ring_prec:
            return coeffs, None
        series_prec = ring_prec

    prec = min(series_prec, ring_prec)
    if deg < prec:
        return coeffs, series_prec

    coeffs = dup_truncate(coeffs, prec, dom)
    return coeffs, prec


def _useries_valuation(s: USeries[Er], dom: Domain) -> int:
    """
    Returns the valuation of this power series.

    If there are no known nonzero coefficients, returns -1.
    """
    coeffs, _ = s

    if not coeffs:
        return -1

    i = -1
    while dom.is_zero(coeffs[i]):
        i -= 1

    return -i - 1


def _unify_prec(
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
) -> tuple[dup[Er], dup[Er], int]:
    """Unify the precision of two series."""
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if prec1 == prec2 == ring_prec:
        return coeffs1, coeffs2, ring_prec
    elif prec1 is not None and prec2 is not None and prec1 == prec2 <= ring_prec:
        return coeffs1, coeffs2, prec1
    elif prec1 == prec2:
        unified_prec = prec1
    elif prec1 is None:
        unified_prec = prec2
    elif prec2 is None:
        unified_prec = prec1
    else:
        unified_prec = min(prec1, prec2)

    if unified_prec is None:
        unified_prec = ring_prec
    else:
        unified_prec = min(unified_prec, ring_prec)

    d1 = dup_degree(coeffs1)
    d2 = dup_degree(coeffs2)

    if d1 >= unified_prec:
        coeffs1 = dup_truncate(coeffs1, unified_prec, dom)
    if d2 >= unified_prec:
        coeffs2 = dup_truncate(coeffs2, unified_prec, dom)

    return coeffs1, coeffs2, unified_prec


def _useries_equality(
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
) -> bool | None:
    """Check if two power series are equal."""
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if prec1 is None and prec2 is None:
        return s1 == s2

    coeffs1, coeffs2, _ = _unify_prec(s1, s2, dom, ring_prec)
    if coeffs1 != coeffs2:
        return False
    return None


def _useries_equal_repr(s1: USeries[Er], s2: USeries[Er]) -> bool:
    return s1 == s2


def _useries_pos(s: USeries[Er], dom: Domain[Er], ring_prec: int) -> USeries[Er]:
    """Return the positive of a power series (which is the same as the series itself)."""
    coeffs, prec = s
    deg = dup_degree(coeffs)

    if prec is None:
        if deg >= ring_prec:
            coeffs = dup_truncate(coeffs, ring_prec, dom)
            prec = ring_prec
    else:
        prec = min(prec, ring_prec)
        if deg >= prec:
            coeffs = dup_truncate(coeffs, prec, dom)

    return coeffs, prec


def _useries_neg(s: USeries[Er], dom: Domain[Er], ring_prec: int) -> USeries[Er]:
    coeffs, prec = s
    neg_coeffs = dup_neg(coeffs, dom)
    s = neg_coeffs, prec
    return _useries_pos(s, dom, ring_prec)


def _useries_add(
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2
    max_degree = max(dup_degree(coeffs1), dup_degree(coeffs2))

    if prec1 is None and prec2 is None and (max_degree < ring_prec):
        series = dup_add(coeffs1, coeffs2, dom)
        return series, None

    coeffs1, coeffs2, min_prec = _unify_prec(s1, s2, dom, ring_prec)
    return dup_add(coeffs1, coeffs2, dom), min_prec


def _useries_add_ground(
    s: USeries[Er], n: Er, dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    """
    Helper function to add a ground element to a power series.

    Mostly need for the series expansion of functions like exp, log, etc.
    """
    coeffs, prec = s

    if n == 0:
        return s

    coeffs = dup_add_ground(coeffs, n, dom)
    return _useries(coeffs, prec, dom, ring_prec)


def _useries_sub(
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2
    max_degree = max(dup_degree(coeffs1), dup_degree(coeffs2))

    if prec1 is None and prec2 is None and (max_degree < ring_prec):
        series = dup_sub(coeffs1, coeffs2, dom)
        return series, None

    coeffs1, coeffs2, min_prec = _unify_prec(s1, s2, dom, ring_prec)
    return dup_sub(coeffs1, coeffs2, dom), min_prec


def _useries_sub_ground(
    s: USeries[Er], n: Er, dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    coeffs, prec = s

    if n == 0:
        return s

    coeffs = dup_sub_ground(coeffs, n, dom)
    return _useries(coeffs, prec, dom, ring_prec)


def _useries_rsub_ground(
    s: USeries[Er], n: Er, dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    """
    Helper function to subtract a power series from the ground element.

    Mostly need for the series expansion of functions like exp, log, etc.
    """
    coeffs, prec = s

    if n == 0:
        return s

    coeffs_neg = dup_neg(coeffs, dom)
    coeffs = dup_add_ground(coeffs_neg, n, dom)
    return _useries(coeffs, prec, dom, ring_prec)


def _useries_mul(
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if prec1 is None and prec2 is None:
        d = dup_degree(coeffs1) + dup_degree(coeffs2)
        if d < ring_prec:
            return dup_mul(coeffs1, coeffs2, dom), None

    coeffs1, coeffs2, min_prec = _unify_prec(s1, s2, dom, ring_prec)
    coeffs = dup_series_mul(coeffs1, coeffs2, min_prec, dom)
    return coeffs, min_prec


def _useries_mul_ground(
    s: USeries[Er], n: Er, dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    coeffs, prec = s

    if n == 0:
        return [], prec
    if n == 1:
        return s
    if n == -1:
        return _useries_neg(s, dom, ring_prec)

    coeffs = dup_mul_ground(coeffs, n, dom)
    return _useries(coeffs, prec, dom, ring_prec)


def _useries_square(s: USeries[Er], dom: Domain[Er], ring_prec: int) -> USeries[Er]:
    """Compute the square of a power series."""
    coeffs, prec = s

    if prec is None:
        d = dup_degree(coeffs) * 2
        if d < ring_prec:
            return dup_series_sqr(coeffs, ring_prec, dom), None
        else:
            prec = ring_prec

    return dup_series_sqr(coeffs, prec, dom), prec


def _useries_sqrt_newton(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the square root of a power series using Newton Iteration."""
    coeffs, prec = s

    if not coeffs:
        return s

    if not dom.is_one(coeffs[-1]):
        raise ValueError("Square root requires the constant term to be one.")

    if prec is None:
        prec = ring_prec

    p: USeries[Ef] = ([dom.one], prec)
    # P_{n+1} = 1/2 * (P_n + S / P_n)
    for precx in _giant_steps(prec):
        p = p[0], precx
        p_inv = _useries_inverse(p, dom, precx)
        tmp = _useries_mul(s, p_inv, dom, precx)
        tmp = _useries_add(p, tmp, dom, precx)
        p = _useries_div_ground(tmp, dom(2), dom, precx)

    return p


def _useries_sqrt(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """
    Computes the square root of a power series.
    """
    coeffs, prec = s

    if not coeffs:
        return s

    c0 = coeffs[-1]

    if dom.is_zero(c0):
        raise ValueError(
            "Cannot compute square root of a series with zero constant term."
        )

    sqrt_c0 = dom.exsqrt(c0)
    if not sqrt_c0:
        raise ValueError("Constant term is not a perfect square in the domain.")

    # Normalize
    norm_s = _useries_div_ground(s, c0, dom, ring_prec)

    if prec is None:
        ds = dup_degree(coeffs)
        if ds % 2 == 0:
            dp = ds // 2
            required_prec = dp + 1
            root = _useries_sqrt_newton(norm_s, dom, required_prec)
            tmp_root = (root[0], ds + 1)
            squared = _useries_square(tmp_root, dom, ds + 1)

            if squared[0] == norm_s[0]:
                root = (root[0], None)
                return _useries_mul_ground(root, sqrt_c0, dom, ring_prec)

    approx_prec = min(prec, ring_prec) if prec is not None else ring_prec
    sqrt_norm_s = _useries_sqrt_newton(norm_s, dom, approx_prec)
    # Denormalize
    sqrt_s = _useries_mul_ground(sqrt_norm_s, sqrt_c0, dom, approx_prec)
    return sqrt_s[0], approx_prec


def _useries_div_direct(
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    """Direct division when divisor has zero valuation."""
    if not dom.is_unit(s2[0][-1]):
        raise ValueError("Trailing coefficient of the divisor must be a unit")

    _, _, min_prec = _unify_prec(s1, s2, dom, ring_prec)
    s2_inv = _useries_inverse(s2, dom, min_prec)
    return _useries_mul(s1, s2_inv, dom, ring_prec)


def _useries_div(
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if not coeffs2:
        raise ZeroDivisionError("Series division by zero")

    val1 = _useries_valuation(s1, dom)
    val2 = _useries_valuation(s2, dom)

    if val1 < val2:
        raise ValueError("quotient would not be a power series")

    if prec1 is None and prec2 is None:
        try:
            q = dup_exquo(coeffs1, coeffs2, dom)
            return q, None
        except ExactQuotientFailed:
            pass

    if val2 == 0:
        return _useries_div_direct(s1, s2, dom, ring_prec)
    else:
        # Shift both series to make divisor's valuation zero
        shifted_coeffs1 = dup_rshift(coeffs1, val2, dom)
        shifted_coeffs2 = dup_rshift(coeffs2, val2, dom)

        prec1 = prec1 - val2 if prec1 is not None else ring_prec
        prec2 = prec2 - val2 if prec2 is not None else ring_prec

        shifted_s1 = _useries(shifted_coeffs1, prec1, dom, ring_prec)
        shifted_s2 = _useries(shifted_coeffs2, prec2, dom, ring_prec)

        return _useries_div_direct(shifted_s1, shifted_s2, dom, ring_prec)


def _useries_div_ground(
    s: USeries[Ef], n: Ef, dom: Domain[Ef], ring_prec: int
) -> USeries[Ef]:
    coeffs, prec = s

    if n == 0:
        raise ZeroDivisionError("Division by zero in power series")
    if n == 1:
        return s
    if n == -1:
        return _useries_neg(s, dom, ring_prec)

    series = dup_mul_ground(coeffs, dom.revert(n), dom)
    return _useries(series, prec, dom, ring_prec)


def _useries_pow_int(
    s: USeries[Er], n: int, dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    """Raise a power series to a integer power with truncation."""
    if n < 0:
        n = -n
        s = _useries_pow_int(s, n, dom, ring_prec)
        try:
            inv = _useries_inverse(s, dom, ring_prec)
            return inv
        except NotReversible:
            raise ValueError("Result would not be a power series")

    coeffs, prec = s

    if n == 0:
        return [dom.one], prec

    if prec is None:
        deg = dup_degree(coeffs) * n
        if deg < ring_prec:
            return dup_series_pow(coeffs, n, ring_prec, dom), None
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    series = dup_series_pow(coeffs, n, prec, dom)
    return series, prec


def _useries_truncate(s: USeries[Er], n: int, dom: Domain[Er]) -> USeries[Er]:
    """Truncate a power series to the first n terms."""
    coeffs, prec = s

    if n < 0:
        raise ValueError("Truncation precision must be non-negative")

    deg = dup_degree(coeffs)
    if deg < n:
        return coeffs, prec
    return dup_truncate(coeffs, n, dom), n


def _useries_compose(
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    """Compose two power series."""
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if coeffs2 and not dom.is_zero(coeffs2[-1]):
        raise ValueError(
            "Series composition requires the constant term of the second series to be zero."
        )

    if prec1 is None and prec2 is None:
        comp = dup_series_compose(coeffs1, coeffs2, ring_prec, dom)

        deg1 = dup_degree(coeffs1)
        deg2 = dup_degree(coeffs2)

        if deg1 * deg2 < ring_prec:
            return comp, None
        else:
            return comp, ring_prec

    coeffs1, coeffs2, min_prec = _unify_prec(s1, s2, dom, ring_prec)
    comp = dup_series_compose(coeffs1, coeffs2, min_prec, dom)
    return comp, min_prec


def _useries_inverse(s: USeries[Er], dom: Domain[Er], ring_prec: int) -> USeries[Er]:
    """Compute the series multiplicative inverse of a power series."""
    coeffs, prec = s

    if not coeffs or not dom.is_unit(coeffs[-1]):
        raise NotReversible("Series inverse requires the constant term to be a unit")

    if prec is None:
        if len(coeffs) == 1:
            return [dom.revert(coeffs[0])], None
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    inv = dup_revert(coeffs, prec, dom)
    inv = dup_truncate(inv, prec, dom)
    return inv, prec


def _useries_reversion(s: USeries[Er], dom: Domain[Er], ring_prec: int) -> USeries[Er]:
    """Compute the composite inverse of a power series."""
    coeffs, prec = s

    if not coeffs or not dom.is_zero(coeffs[-1]):
        raise NotReversible(
            "Series compositional inverse requires the constant term to be zero."
        )

    if len(coeffs) >= 2 and not dom.is_unit(coeffs[-2]):
        raise NotReversible(
            "Series compositional inverse requires the linear term to be unit."
        )

    if prec is None:
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    series = dup_series_reversion(coeffs, prec, dom)
    return series, prec


def _useries_derivative(s: USeries[Er], dom: Domain[Er], ring_prec: int) -> USeries[Er]:
    """Compute the first derivative of a power series."""
    coeffs, prec = s
    series = dup_diff(coeffs, 1, dom)
    if prec:
        prec -= 1
    return _useries(series, prec, dom, ring_prec)


def _useries_integrate(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the integral of a power series."""
    coeffs, prec = s
    series = dup_integrate(coeffs, 1, dom)
    if prec:
        prec += 1
    return _useries(series, prec, dom, ring_prec)


def _useries_log(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the logarithm of a power series."""
    coeffs, prec = s

    if not coeffs:
        raise ValueError("Logarithm of a zero series is undefined.")

    if not dom.is_one(coeffs[-1]):
        raise ValueError(
            "Logarithm of a power series requires the constant term to be one."
        )

    if len(coeffs) == 1:
        return ([], prec)

    if prec is None:
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    s = coeffs, prec

    if len(coeffs) > 20:
        # log(f(x))' = f'(x) / f(x)
        dv_s = _useries_derivative(s, dom, ring_prec)
        inverse_s = _useries_inverse(s, dom, ring_prec)
        log_derivative = _useries_mul(dv_s, inverse_s, dom, ring_prec)
        return _useries_integrate(log_derivative, dom, ring_prec)

    else:
        c: list[Ef] = [dom.zero]
        for k in range(1, prec):
            c.append(dom.quo((-1) ** (k - 1), dom(k)))

        f: USeries[Ef] = c[::-1], prec
        s = _useries_add_ground(s, dom(-1), dom, ring_prec)
        r = _useries_compose(f, s, dom, ring_prec)
        return r


def _useries_log1p(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Computes log(1 + s) of power series with zero constant term."""

    if not dom.is_zero(s[0][-1]):
        raise ValueError(
            "Log1p requires the constant term of the input series to be zero."
        )

    s = _useries_add_ground(s, dom.one, dom, ring_prec)
    return _useries_log(s, dom, ring_prec)


def _useries_exp(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    coeffs, prec = s

    if not coeffs:
        return ([dom.one], prec)

    if not dom.is_zero(coeffs[-1]):
        raise ValueError(
            "Exponential requires the constant term of the input series to be zero."
        )

    if prec is None:
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    s = coeffs, prec

    one: USeries[Ef] = ([dom.one], prec)

    if len(coeffs) > 20:
        exp = one

        for precx in _giant_steps(prec):
            exp = exp[0], precx
            log_exp = _useries_log(exp, dom, ring_prec)
            diff = _useries_sub(s, log_exp, dom, ring_prec)
            corr = _useries_mul(diff, exp, dom, ring_prec)
            exp = _useries_add(exp, corr, dom, ring_prec)

        return exp

    else:
        n: Ef = dom(1)
        c: list[Ef] = []
        for k in range(prec):
            c.append(dom.revert(n))
            k += 1
            n *= dom(k)

        f: USeries[Ef] = c[::-1], prec
        r = _useries_compose(f, s, dom, ring_prec)
        return r


def _useries_expm1(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Computes exp(s) - 1 of power series."""
    exp_s = _useries_exp(s, dom, ring_prec)
    return _useries_add_ground(exp_s, dom(-1), dom, ring_prec)


def _useries_atan(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the arctangent of a power series."""
    coeffs, prec = s

    if not coeffs:
        return s

    if not dom.is_zero(coeffs[-1]):
        raise ValueError(
            "Arctangent requires the constant term of the input series to be zero."
        )

    if prec is None:
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    s = coeffs, prec

    # atan(f(x))' = f'(x) / (1 + f(x)**2)
    ds = _useries_derivative(s, dom, ring_prec)
    s2 = _useries_square(s, dom, ring_prec)
    s2 = _useries_add_ground(s2, dom.one, dom, ring_prec)  # 1 + s^2
    inv_s2 = _useries_inverse(s2, dom, ring_prec)
    dv_s = _useries_mul(ds, inv_s2, dom, ring_prec)
    return _useries_integrate(dv_s, dom, ring_prec)


def _useries_atanh(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the hyperbolic arctanh of a power series"""
    coeffs, prec = s

    if not coeffs:
        return s

    if not dom.is_zero(coeffs[-1]):
        raise ValueError(
            "Hyperbolic Arctangent requires the constant term of the input series to be zero."
        )

    if prec is None:
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    s = coeffs, prec

    # atanh(f(x))' = f'(x) / (1 - f(x)**2)
    ds = _useries_derivative(s, dom, ring_prec)
    s2 = _useries_square(s, dom, ring_prec)
    q = _useries_rsub_ground(s2, dom.one, dom, ring_prec)
    q_inv = _useries_inverse(q, dom, ring_prec)
    dv_atanh = _useries_mul(ds, q_inv, dom, ring_prec)
    return _useries_integrate(dv_atanh, dom, ring_prec)


def _useries_asin(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the arcsine of a power series."""
    coeffs, prec = s

    if not coeffs:
        return s

    if not dom.is_zero(coeffs[-1]):
        raise ValueError(
            "Hyperbolic arcsine requires the constant term of the input series to be zero."
        )

    if prec is None:
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    s = coeffs, prec

    if len(coeffs) > 20:
        # asin(f(x))' = f'(x) / sqrt(1 - f(x)**2)
        ds = _useries_derivative(s, dom, ring_prec)
        s2 = _useries_square(s, dom, ring_prec)
        q = _useries_rsub_ground(s2, dom.one, dom, ring_prec)
        q_sqrt = _useries_sqrt(q, dom, ring_prec)
        q_sqrt_inv = _useries_inverse(q_sqrt, dom, ring_prec)  # 1 / sqrt(1 - s^2)

        dv_asin = _useries_mul(ds, q_sqrt_inv, dom, ring_prec)
        return _useries_integrate(dv_asin, dom, ring_prec)

    else:
        c: list[Ef] = [dom.zero, dom.one, dom.zero]
        for k in range(3, prec, 2):
            term = dom(k - 2) ** 2 * c[-2]
            term = dom.quo(term, dom(k * (k - 1)))
            c.append(term)
            c.append(dom.zero)

        f: USeries[Ef] = (c[::-1], prec)
        return _useries_compose(f, s, dom, ring_prec)


def _useries_asinh(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the hyperbolic arcsin of a power series"""
    coeffs, prec = s

    if not coeffs:
        return s

    if not dom.is_zero(coeffs[-1]):
        raise ValueError(
            "Hyperbolic Arcsine requires the constant term of the input series to be zero."
        )

    if prec is None:
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    s = coeffs, prec

    if len(coeffs) > 20:
        # asinh(f(x))' = f'(x) / sqrt(1 + f(x)**2)
        ds = _useries_derivative(s, dom, ring_prec)
        s2 = _useries_square(s, dom, ring_prec)
        q = _useries_add_ground(s2, dom.one, dom, ring_prec)
        q_sqrt = _useries_sqrt(q, dom, ring_prec)
        q_sqrt_inv = _useries_inverse(q_sqrt, dom, ring_prec)  # 1 / sqrt(1 + s^2)

        dv_asinh = _useries_mul(ds, q_sqrt_inv, dom, ring_prec)
        return _useries_integrate(dv_asinh, dom, ring_prec)

    else:
        c: list[Ef] = [dom.zero, dom.one, dom.zero]

        for k in range(3, prec, 2):
            term = (k - 2) ** 2 * c[-2]
            term = dom.quo(-term, dom(k * (k - 1)))
            c.append(term)
            c.append(dom.zero)

        f: USeries[Ef] = c[::-1], prec
        return _useries_compose(f, s, dom, ring_prec)


def _useries_tan(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the tangent of a power series."""
    coeffs, prec = s

    if not coeffs:
        return s

    if not dom.is_zero(coeffs[-1]):
        raise ValueError(
            "Tangent requires the constant term of the input series to be zero."
        )

    if prec is None:
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    s = coeffs, prec

    y: USeries[Ef] = ([], prec)
    # y_{n+1} = y_n + (s - \arctan(y_n)) * (1 + y_n^2)
    #                       term1            term2
    for precx in _giant_steps(prec):
        y = y[0], precx
        term1 = _useries_atan(y, dom, ring_prec)
        term1 = _useries_sub(s, term1, dom, ring_prec)

        term2 = _useries_square(y, dom, ring_prec)
        term2 = _useries_add_ground(term2, dom.one, dom, ring_prec)

        mul = _useries_mul(term1, term2, dom, ring_prec)
        y = _useries_add(y, mul, dom, ring_prec)

    return y


def _useries_tanh(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the hyperbolic tangent of a power series."""
    coeffs, prec = s

    if not coeffs:
        return s

    if not dom.is_zero(coeffs[-1]):
        raise ValueError(
            "Hyperbolic Tangent requires the constant term of the input series to be zero."
        )

    if prec is None:
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    s = coeffs, prec

    y: USeries[Ef] = ([], prec)
    # y_{n+1} = y_n + (s - \arctanh(y_n)) * (1 - y_n^2)
    #                       term1            term2
    for precx in _giant_steps(prec):
        y = y[0], precx
        term1 = _useries_atanh(y, dom, ring_prec)
        term1 = _useries_sub(s, term1, dom, ring_prec)

        term2 = _useries_square(y, dom, ring_prec)
        term2 = _useries_rsub_ground(term2, dom.one, dom, ring_prec)

        mul = _useries_mul(term1, term2, dom, ring_prec)
        y = _useries_add(y, mul, dom, ring_prec)

    return y


def _useries_sin(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the sine of a power series."""
    coeffs, prec = s

    if not coeffs:
        return s

    if not dom.is_zero(coeffs[-1]):
        raise ValueError(
            "Sine requires the constant term of the input series to be zero."
        )

    if prec is None:
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    s = coeffs, prec

    if len(coeffs) > 20:
        # sin(f(x)) = 2 * tan(f(x)/2) / (1 + tan(f(x)/2)**2)
        s = _useries_div_ground(s, dom(2), dom, ring_prec)
        # t = tan(s/2)
        t = _useries_tan(s, dom, ring_prec)
        t_sq = _useries_square(t, dom, ring_prec)

        p = _useries_mul_ground(t, dom(2), dom, ring_prec)  # 2*t

        q = _useries_add_ground(t_sq, dom.one, dom, ring_prec)
        q_inv = _useries_inverse(q, dom, ring_prec)  # 1/(1 + t^2)

        return _useries_mul(p, q_inv, dom, ring_prec)

    else:
        n: Ef = dom.one
        c: list[Ef] = [dom.zero]

        for k in range(2, prec + 2, 2):
            c.append(dom.revert(n))
            c.append(dom.zero)
            n *= -k * (k + 1)

        f: USeries[Ef] = (c[::-1], prec)
        return _useries_compose(f, s, dom, ring_prec)


def _useries_sinh(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the hyperbolic sine of a power series"""
    coeffs, prec = s

    if not coeffs:
        return s

    if not dom.is_zero(coeffs[-1]):
        raise ValueError(
            "Hyperbolic sine requires the constant term of the input series to be zero."
        )

    if prec is None:
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    s = coeffs, prec

    if len(coeffs) > 40:
        # sinh(f(x)) = (exp(f(x)) - exp(-f(x))) / 2
        e = _useries_exp(s, dom, ring_prec)
        e_inv = _useries_inverse(e, dom, ring_prec)
        diff = _useries_sub(e, e_inv, dom, ring_prec)
        return _useries_div_ground(diff, dom(2), dom, ring_prec)

    else:
        n: Ef = dom.one
        c: list[Ef] = [dom.zero]

        for k in range(2, prec + 2, 2):
            c.append(dom.revert(n))
            c.append(dom.zero)
            n *= k * (k + 1)

        f: USeries[Ef] = c[::-1], prec

        return _useries_compose(f, s, dom, ring_prec)


def _useries_cos(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the cosine of a power series."""
    coeffs, prec = s

    if not coeffs:
        return ([dom.one], prec)

    if not dom.is_zero(coeffs[-1]):
        raise ValueError(
            "Cosine requires the constant term of the input series to be zero."
        )

    if prec is None:
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    s = coeffs, prec

    if len(coeffs) > 20:
        # cos(f(x)) = (1 - tan(f(x)/2)**2) / (1 + tan(f(x)/2)**2)
        s = _useries_div_ground(s, dom(2), dom, ring_prec)
        # t = tan(s/2)
        t = _useries_tan(s, dom, ring_prec)

        t_sq = _useries_square(t, dom, ring_prec)  # t_sq = tan^2(s/2)
        p = _useries_rsub_ground(t_sq, dom.one, dom, ring_prec)  # 1 - tan^2(s/2)

        q = _useries_add_ground(t_sq, dom.one, dom, ring_prec)
        q_inv = _useries_inverse(q, dom, ring_prec)  # 1/(1 + tan^2(s/2))

        return _useries_mul(p, q_inv, dom, ring_prec)
    else:
        n: Ef = dom.one
        c: list[Ef] = []

        for k in range(2, prec + 2, 2):
            c.append(dom.revert(n))
            c.append(dom.zero)
            n *= -k * (k - 1)

        f: USeries[Ef] = (c[::-1], prec)
        return _useries_compose(f, s, dom, ring_prec)


def _useries_cosh(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the hyperbolic cosine of a power series"""
    coeffs, prec = s

    if not coeffs:
        return ([dom.one], prec)

    if not dom.is_zero(coeffs[-1]):
        raise ValueError(
            "Hyperbolic cosine requires the constant term of the input series to be zero."
        )

    if prec is None:
        prec = ring_prec
    else:
        prec = min(prec, ring_prec)

    s = coeffs, prec

    if len(coeffs) > 40:
        # cosh(f(x)) = (exp(f(x)) + exp(-f(x))) / 2
        e = _useries_exp(s, dom, ring_prec)
        e_inv = _useries_inverse(e, dom, ring_prec)
        diff = _useries_add(e, e_inv, dom, ring_prec)
        return _useries_div_ground(diff, dom(2), dom, ring_prec)
    else:
        n: Ef = dom.one
        c: list[Ef] = [dom.one, dom.zero]

        for k in range(2, prec, 2):
            n *= (k - 1) * k
            c.append(dom.revert(n))
            c.append(dom.zero)

        f: USeries[Ef] = c[::-1], prec

        return _useries_compose(f, s, dom, ring_prec)


def _useries_hypot(
    s1: USeries[Ef], s2: USeries[Ef], dom: Field[Ef], ring_prec: int
) -> USeries[Ef]:
    """Compute the hypotenuse of two power series."""
    # hypot(f(x), g(x)) = sqrt(f(x)**2 + g(x)**2)
    sqr_s1 = _useries_square(s1, dom, ring_prec)
    sqr_s2 = _useries_square(s2, dom, ring_prec)
    add = _useries_add(sqr_s1, sqr_s2, dom, ring_prec)
    return _useries_sqrt(add, dom, ring_prec)


class PythonPowerSeriesRingZZ:
    """
    Python implementation of power series ring over integers :ref:`ZZ`.

    This class provides comprehensive power series operations over the integer ring,
    supporting both series manipulations with precision handling and truncation.

    Parameters
    ==========

    prec : int, optional
        The default precision for power series operations. Default is 6.

    Examples
    ========

    >>> from sympy.polys.series.ringpython import PythonPowerSeriesRingZZ
    >>> R = PythonPowerSeriesRingZZ()
    >>> s = R([1, 2, 3])  # 1 + 2*x + 3*x^2
    >>> R.print(s)
    1 + 2*x + 3*x**2

    >>> s_pow = R.pow_int(s, 2)  # Square the series
    >>> R.print(s_pow)
    1 + 4*x + 10*x**2 + 12*x**3 + 9*x**4

    >>> s_inv = R.inverse(R([1, 1]))  # Inverse of 1 + x
    >>> R.print(s_inv)
    1 - x + x**2 - x**3 + x**4 - x**5 + O(x**6)

    Note
    ====

    The recommended way to create a power series ring is using the factory function
    which returns a new instance of the higher level PowerSeriesRing class with
    the ring generator:

    >>> from sympy.polys.series import power_series_ring
    >>> from sympy import ZZ
    >>> R, x = power_series_ring("x", ZZ, 6)
    >>> R
    Power Series Ring in x over ZZ of size 6
    >>> type(x)
    <class 'sympy.polys.series.ring.PowerSeriesElement'>

    This function automatically uses the Flint implementation if available for better
    performance, falling back to the Python implementation otherwise.

    See Also
    ========

    sympy.polys.series.ringpython.PythonPowerSeriesRingQQ
    sympy.polys.series.ringflint.FlintPowerSeriesRingZZ
    sympy.polys.series.ring.power_series_ring
    sympy.polys.series.ring.PowerSeriesRingRing
    sympy.polys.series.ring.PowerSeriesRingField
    sympy.polys.series.ring.PowerSeriesElement
    """

    _domain = ZZ

    def __init__(self, prec: int = 6) -> None:
        if prec < 0:
            raise ValueError("Power series precision must be non-negative")
        self._prec = prec

    def __repr__(self) -> str:
        return (
            f"Python Power Series Ring over {self._domain} with precision {self._prec}"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, PythonPowerSeriesRingZZ):
            return NotImplemented
        return self._prec == other.prec

    def __hash__(self) -> int:
        return hash((self._domain, self._prec))

    def __call__(
        self, coeffs: Sequence[MPZ | int], prec: int | None = None
    ) -> USeries[MPZ]:
        """
        Create a power series from a list of coefficients.

        If `prec` is not specified, it defaults to the ring's precision.
        """
        s: list[MPZ] = []
        for c in coeffs:
            if isinstance(c, MPZ):
                s.append(c)
            elif isinstance(c, int):
                s.append(self._domain(c))
            else:
                raise TypeError(f"Unsupported coefficient type: {type(c)}")

        return self.from_list(s, prec)

    @property
    def domain(self) -> Domain[MPZ]:
        """Return the ground domain of the power series ring."""
        return self._domain

    @property
    def prec(self) -> int:
        """Return the ring's precision."""
        return self._prec

    @property
    def one(self) -> USeries[MPZ]:
        if self._prec == 0:
            return ([], 0)
        return ([self._domain.one], None)

    @property
    def zero(self) -> USeries[MPZ]:
        if self._prec == 0:
            return ([], 0)
        return ([], None)

    @property
    def gen(self) -> USeries[MPZ]:
        if self._prec < 2:
            return ([], self._prec)
        return ([self._domain.one, self._domain.zero], None)

    def pretty(
        self, s: USeries[MPZ], *, symbol: str = "x", ascending: bool = True
    ) -> str:
        """Return a pretty-printed string representation of a power series."""
        coeffs, prec = s
        return series_pprint(coeffs, prec, sym=symbol, ascending=ascending)

    def print(
        self, s: USeries[MPZ], *, symbol: str = "x", ascending: bool = True
    ) -> None:
        """Print a pretty-printed representation of a power series."""
        print(self.pretty(s, symbol=symbol, ascending=ascending))

    def from_list(self, coeffs: list[MPZ], prec: int | None = None) -> USeries[MPZ]:
        """
        Create a power series from a list of ground coefficients.

        If `prec` is not specified, it defaults to the ring's precision.
        """
        coeffs = dup_reverse(coeffs, ZZ)
        if prec is None:
            if len(coeffs) <= self._prec:
                return coeffs, None
            else:
                prec = self._prec
        else:
            prec = min(prec, self._prec)

        if len(coeffs) > prec:
            coeffs = dup_truncate(coeffs, prec, self._domain)

        return coeffs, prec

    def to_list(self, s: USeries[MPZ]) -> list[MPZ]:
        """Returns the list of series coefficients."""
        coeffs, _ = s
        return coeffs[::-1]

    def to_dense(self, s: USeries[MPZ]) -> dup[MPZ]:
        """Return the coefficients of a power series as a dense list."""
        return list(s[0])

    def series_prec(self, s: USeries[MPZ]) -> int | None:
        """Return the precision of a power series."""
        _, prec = s
        return prec

    def equal(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> bool | None:
        """Check if two power series are equal up to their minimum precision."""
        return _useries_equality(s1, s2, self._domain, self._prec)

    def equal_repr(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> bool:
        """Check if two power series have the same representation."""
        return _useries_equal_repr(s1, s2)

    def is_ground(self, arg: USeries[MPZ]) -> bool | None:
        """Check if a arg is a ground element of the power series ring."""
        if self.prec == 0:
            return None

        return len(self.to_list(arg)) <= 1

    def constant_coefficient(self, s: USeries[MPZ]) -> MPZ:
        """Return the constant coefficient of a power series."""
        coeffs, _ = s
        if len(coeffs) > 0:
            return coeffs[-1]
        return self._domain.zero

    def positive(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """Return the unary positive of a power series, adjusted to the ring's precision."""
        return _useries_pos(s, self._domain, self._prec)

    def negative(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """Return the unary negative of a power series."""
        return _useries_neg(s, self._domain, self._prec)

    def add(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """Add two power series."""
        return _useries_add(s1, s2, self._domain, self._prec)

    def add_ground(self, s: USeries[MPZ], n: MPZ) -> USeries[MPZ]:
        """Add a ground element to a power series."""
        return _useries_add_ground(s, n, self._domain, self._prec)

    def subtract(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """Subtract two power series."""
        return _useries_sub(s1, s2, self._domain, self._prec)

    def subtract_ground(self, s: USeries[MPZ], n: MPZ) -> USeries[MPZ]:
        """Subtract a ground element from a power series."""
        return _useries_sub_ground(s, n, self._domain, self._prec)

    def rsubtract_ground(self, s: USeries[MPZ], n: MPZ) -> USeries[MPZ]:
        """Subtract a power series from a ground element."""
        return _useries_rsub_ground(s, n, self._domain, self._prec)

    def multiply(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """Multiply two power series."""
        return _useries_mul(s1, s2, self._domain, self._prec)

    def multiply_ground(self, s: USeries[MPZ], n: MPZ) -> USeries[MPZ]:
        """Multiply a power series by a ground element."""
        return _useries_mul_ground(s, n, self._domain, self._prec)

    def divide(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """Divide two power series."""
        return _useries_div(s1, s2, self._domain, self._prec)

    def pow_int(self, s: USeries[MPZ], n: int) -> USeries[MPZ]:
        """Raise a power series to a integer power."""
        return _useries_pow_int(s, n, self._domain, self._prec)

    def square(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """Compute the square of a power series."""
        return _useries_mul(s, s, self._domain, self._prec)

    def compose(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """Compose two power series, `s1(s2)`."""
        return _useries_compose(s1, s2, self._domain, self._prec)

    def inverse(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """Compute the multiplicative inverse of a power series."""
        return _useries_inverse(s, self._domain, self._prec)

    def reversion(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """Compute the compositional inverse of a power series."""
        return _useries_reversion(s, self._domain, self._prec)

    def truncate(self, s: USeries[MPZ], n: int) -> USeries[MPZ]:
        """Truncate a power series to `n` terms."""
        return _useries_truncate(s, n, self._domain)

    def differentiate(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """Compute the derivative of a power series."""
        return _useries_derivative(s, self._domain, self._prec)


class PythonPowerSeriesRingQQ:
    """
    Python implementation of power series ring over rational field :ref:`QQ`.

    This class provides comprehensive power series operations over the rational field,
    supporting series manipulations with precision handling and truncation.
    It extends the integer ring functionality with support for rational coefficients
    and integration.

    Parameters
    ==========

    prec : int, optional
        The default precision for power series operations. Default is 6.

    Examples
    ========

    >>> from sympy.polys.series.ringpython import PythonPowerSeriesRingQQ
    >>> R = PythonPowerSeriesRingQQ()
    >>> s = R([1, (1, 2), (1, 3)])  # 1 + x/2 + x^2/3
    >>> R.print(s)
    1 + 1/2*x + 1/3*x**2

    >>> s_int = R.integrate(s)  # Integration
    >>> R.print(s_int)
    x + 1/4*x**2 + 1/9*x**3

    >>> s_inv = R.inverse(R([1, (1, 2)]))  # Inverse of 1 + x/2
    >>> R.print(s_inv)
    1 - 1/2*x + 1/4*x**2 - 1/8*x**3 + 1/16*x**4 - 1/32*x**5 + O(x**6)

    Note
    ====

    The recommended way to create a power series ring is using the factory function
    which returns a new instance of the higher level PowerSeriesRing class with
    the ring generator:

    >>> from sympy.polys.series import power_series_ring
    >>> from sympy import QQ
    >>> R, x = power_series_ring("x", QQ, 6)
    >>> R
    Power Series Ring in x over QQ of size 6
    >>> type(x)
    <class 'sympy.polys.series.ring.PowerSeriesElement'>

    This function automatically uses the Flint implementation if available for better
    performance, falling back to the Python implementation otherwise.

    See Also
    ========

    sympy.polys.series.ringpython.PythonPowerSeriesRingZZ
    sympy.polys.series.ringflint.FlintPowerSeriesRingQQ
    sympy.polys.series.ring.power_series_ring
    sympy.polys.series.ring.PowerSeriesRingRing
    sympy.polys.series.ring.PowerSeriesRingField
    sympy.polys.series.ring.PowerSeriesElement
    """

    _domain = QQ

    def __init__(self, prec: int = 6) -> None:
        if prec < 0:
            raise ValueError("Power series precision must be non-negative")
        self._prec = prec

    def __repr__(self) -> str:
        return (
            f"Python Power Series Ring over {self._domain} with precision {self._prec}"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, PythonPowerSeriesRingQQ):
            return NotImplemented
        return self._prec == other.prec

    def __hash__(self) -> int:
        return hash((self._domain, self._prec))

    def __call__(
        self, coeffs: Sequence[MPQ | int | tuple[int, int]], prec: int | None = None
    ) -> USeries[MPQ]:
        """
        Create a power series from a list of coefficients.

        If `prec` is not specified, it defaults to the ring's precision.
        """
        s: list[MPQ] = []
        for c in coeffs:
            if isinstance(c, MPQ):
                s.append(c)
            elif isinstance(c, int):
                s.append(self._domain(c))
            elif isinstance(c, tuple):
                s.append(self._domain(*c))
            else:
                raise TypeError(f"Unsupported coefficient type: {type(c)}")

        return self.from_list(s, prec)

    @property
    def domain(self) -> Domain[MPQ]:
        """Return the ground domain of the power series ring."""
        return self._domain

    @property
    def prec(self) -> int:
        """Return the ring's precision."""
        return self._prec

    @property
    def one(self) -> USeries[MPQ]:
        if self._prec == 0:
            return ([], 0)
        return ([QQ(1)], None)

    @property
    def zero(self) -> USeries[MPQ]:
        if self._prec == 0:
            return ([], 0)
        return ([], None)

    @property
    def gen(self) -> USeries[MPQ]:
        if self._prec < 2:
            return ([], self._prec)
        return ([self._domain.one, self._domain.zero], None)

    def pretty(
        self, s: USeries[MPQ], *, symbol: str = "x", ascending: bool = True
    ) -> str:
        """Return a pretty-printed string representation of a power series."""
        coeffs, prec = s
        return series_pprint(coeffs, prec, sym=symbol, ascending=ascending)

    def print(
        self, s: USeries[MPQ], *, symbol: str = "x", ascending: bool = True
    ) -> None:
        """Print a pretty-printed representation of a power series."""
        print(self.pretty(s, symbol=symbol, ascending=ascending))

    def from_list(self, coeffs: list[MPQ], prec: int | None = None) -> USeries[MPQ]:
        """
        Create a power series from a list of ground coefficients.

        If `prec` is not specified, it defaults to the ring's precision.
        """
        coeffs = dup_reverse(coeffs, QQ)
        if prec is None:
            if len(coeffs) <= self._prec:
                return coeffs, None
            else:
                prec = self._prec
        else:
            prec = min(prec, self._prec)

        if len(coeffs) > prec:
            coeffs = dup_truncate(coeffs, prec, self._domain)

        return coeffs, prec

    def to_list(self, s: USeries[MPQ]) -> list[MPQ]:
        """Return the list of series coefficients."""
        coeffs, _ = s
        return coeffs[::-1]

    def to_dense(self, s: USeries[MPQ]) -> dup[MPQ]:
        """Return the coefficients of a power series as a dense list."""
        return list(s[0])

    def series_prec(self, s: USeries[MPQ]) -> int | None:
        """Return the precision of a power series."""
        _, prec = s
        return prec

    def equal(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> bool | None:
        """Check if two power series are equal up to their minimum precision."""
        return _useries_equality(s1, s2, self._domain, self._prec)

    def equal_repr(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> bool:
        """Check if two power series have the same representation."""
        return _useries_equal_repr(s1, s2)

    def is_ground(self, arg: USeries[MPQ]) -> bool | None:
        """Check if a arg is a ground element of the power series ring."""
        if self.prec == 0:
            return None

        return len(self.to_list(arg)) <= 1

    def constant_coefficient(self, s: USeries[MPQ]) -> MPQ:
        """Return the constant coefficient of a power series."""
        coeffs, _ = s
        if len(coeffs) > 0:
            return coeffs[-1]
        return self._domain.zero

    def positive(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Return the unary positive of a power series, adjusted to the ring's precision."""
        return _useries_pos(s, self._domain, self._prec)

    def negative(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Return the unary negative of a power series."""
        return _useries_neg(s, self._domain, self._prec)

    def add(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """Add two power series."""
        return _useries_add(s1, s2, self._domain, self._prec)

    def add_ground(self, s: USeries[MPQ], n: MPQ) -> USeries[MPQ]:
        """Add a ground element to a power series."""
        return _useries_add_ground(s, n, self._domain, self._prec)

    def subtract(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """Subtract two power series."""
        return _useries_sub(s1, s2, self._domain, self._prec)

    def subtract_ground(self, s: USeries[MPQ], n: MPQ) -> USeries[MPQ]:
        """Subtract a ground element from a power series."""
        return _useries_sub_ground(s, n, self._domain, self._prec)

    def rsubtract_ground(self, s: USeries[MPQ], n: MPQ) -> USeries[MPQ]:
        """Subtract a power series from a ground element."""
        return _useries_rsub_ground(s, n, self._domain, self._prec)

    def multiply(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """Multiply two power series."""
        return _useries_mul(s1, s2, self._domain, self._prec)

    def multiply_ground(self, s: USeries[MPQ], n: MPQ) -> USeries[MPQ]:
        """Multiply a power series by a ground element."""
        return _useries_mul_ground(s, n, self._domain, self._prec)

    def divide(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """Divide two power series."""
        return _useries_div(s1, s2, self._domain, self._prec)

    def pow_int(self, s: USeries[MPQ], n: int) -> USeries[MPQ]:
        """Raise a power series to a integer power."""
        return _useries_pow_int(s, n, self._domain, self._prec)

    def square(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the square of a power series."""
        return _useries_square(s, self._domain, self._prec)

    def sqrt(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the sqrt of a power series."""
        return _useries_sqrt(s, self._domain, self._prec)

    def compose(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """Compose two power series, `s1(s2)`."""
        return _useries_compose(s1, s2, self._domain, self._prec)

    def inverse(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the multiplicative inverse of a power series."""
        return _useries_inverse(s, self._domain, self._prec)

    def reversion(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the compositional inverse of a power series."""
        return _useries_reversion(s, self._domain, self._prec)

    def truncate(self, s: USeries[MPQ], n: int) -> USeries[MPQ]:
        """Truncate a power series to `n` terms."""
        return _useries_truncate(s, n, self._domain)

    def differentiate(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the derivative of a power series."""
        return _useries_derivative(s, self._domain, self._prec)

    def integrate(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the integral of a power series."""
        return _useries_integrate(s, self._domain, self._prec)

    def log(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the logarithm of a power series."""
        return _useries_log(s, self._domain, self._prec)

    def log1p(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the logarithm of (1 + x) for a power series."""
        return _useries_log1p(s, self._domain, self._prec)

    def exp(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the exponential of a power series."""
        return _useries_exp(s, self._domain, self._prec)

    def expm1(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the exponential of a power series minus 1."""
        return _useries_expm1(s, self._domain, self._prec)

    def atan(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the arctangent of a power series."""
        return _useries_atan(s, self._domain, self._prec)

    def atanh(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the hyperbolic arctangent of a power series."""
        return _useries_atanh(s, self._domain, self._prec)

    def asin(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the arcsine of a power series."""
        return _useries_asin(s, self._domain, self._prec)

    def asinh(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the hyperbolic arcsine of a power series."""
        return _useries_asinh(s, self._domain, self._prec)

    def tan(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the tangent of a power series."""
        return _useries_tan(s, self._domain, self._prec)

    def tanh(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the hyperbolic tangent of a power series."""
        return _useries_tanh(s, self._domain, self._prec)

    def sin(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the sine of a power series."""
        return _useries_sin(s, self._domain, self._prec)

    def sinh(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the hyperbolic sine of a power series."""
        return _useries_sinh(s, self._domain, self._prec)

    def cos(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the cosine of a power series."""
        return _useries_cos(s, self._domain, self._prec)

    def cosh(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the hyperbolic cosine of a power series."""
        return _useries_cosh(s, self._domain, self._prec)

    def hypot(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the hypotenuse of two power series."""
        return _useries_hypot(s1, s2, self._domain, self._prec)
