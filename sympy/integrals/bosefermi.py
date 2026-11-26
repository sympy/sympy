r"""
Definite integrals of Bose-Einstein and Fermi-Dirac type.

This module evaluates integrals of the form

.. math::

    \int_0^\infty \frac{x^p}{e^{a x + b} \pm 1}\,\mathrm{d}x,

which appear throughout statistical physics (Planck's law, the Fermi gas,
white dwarfs, ...).  The result is expressed with the gamma function and the
polylogarithm; for ``b = 0`` the polylogarithm reduces to the Riemann zeta
function (bosons) or to the Dirichlet eta function (fermions).

``meijerint`` cannot compute these integrals: the Mellin transform of
``1/(exp(a*x + b) - 1)`` is ``a**(-s)*gamma(s)*zeta(s)`` (and similarly with
the eta function), which is not a Meijer G function.  Therefore
``_meijerint_definite_2`` calls :func:`bose_fermi_definite` before trying the
Meijer G machinery.
"""

from __future__ import annotations

from sympy.core import Add, S, Wild
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.special.gamma_functions import gamma
from sympy.functions.special.zeta_functions import (
    dirichlet_eta, polylog, zeta)
from sympy.logic.boolalg import And


__all__ = ["bose_fermi_definite"]


def _match_denominator(den, x):
    """Split ``A*exp(a*x + b) + C`` into ``(a, b, sign, scale)``.

    ``sign`` is +1 when ``A/C > 0`` and -1 when ``A/C < 0``, while ``scale``
    is the overall constant ``1/C`` (or ``-1/C``).  ``None`` is returned when
    ``den`` has a different shape; ``A`` and ``C`` must be numbers.
    """
    terms = Add.make_args(den)
    if len(terms) != 2:
        return None

    exp_term = const_term = None
    for term in terms:
        if isinstance(term.as_coeff_Mul()[1], exp):
            if exp_term is not None:
                return None
            exp_term = term
        elif term.is_number:
            if const_term is not None:
                return None
            const_term = term
        else:
            return None

    if exp_term is None or const_term is None:
        return None

    coefficient, exp_part = exp_term.as_coeff_Mul()
    if not coefficient.is_number or coefficient.is_zero:
        return None

    alpha = Wild('alpha', exclude=[x])
    beta = Wild('beta', exclude=[x])
    match = exp_part.args[0].match(alpha*x + beta)
    if match:
        a = match[alpha]
        b = match.get(beta, S.Zero)
    elif exp_part.args[0] == x:
        a, b = S.One, S.Zero
    else:
        return None

    ratio = coefficient/const_term
    if ratio.is_positive:
        return a, b + log(ratio), S.One, S.One/const_term
    if ratio.is_negative:
        return a, b + log(-ratio), S.NegativeOne, -S.One/const_term
    return None


def bose_fermi_definite(f, x):
    """Evaluate ``Integral(f, (x, 0, oo))`` for a Bose-Fermi integrand.

    Returns ``(result, cond)`` in the same way as
    :func:`sympy.integrals.meijerint.meijerint_definite`, or ``None`` when
    ``f`` is not of the supported form or when the integral diverges.

    Examples
    ========

    >>> from sympy import exp, integrate, oo, symbols
    >>> x = symbols('x', positive=True)
    >>> integrate(x**3/(exp(x) - 1), (x, 0, oo))
    pi**4/15
    """
    num, den = f.as_numer_denom()
    coefficient, num = num.as_coeff_Mul()

    if num == 1:
        p = S.Zero
    elif num == x:
        p = S.One
    elif num.is_Pow and num.base == x:
        p = num.exp
    else:
        return None

    den_match = _match_denominator(den, x)
    if den_match is None:
        return None
    a, b, sign, scale = den_match

    cond = And(a > 0, p > -1)
    if sign is S.NegativeOne:
        cond = And(cond, p > 0, b >= 0)
    if cond is S.false:
        return None

    prefactor = coefficient*scale*gamma(p + 1)/a**(p + 1)
    if sign is S.One:
        # fermions: 1/(e^(a*x+b) + 1)
        if b == 0:
            special = dirichlet_eta(p + 1)
        elif p == 0:
            special = log(1 + exp(-b))
        else:
            special = -polylog(p + 1, -exp(-b))
    else:
        # bosons: 1/(e^(a*x+b) - 1)
        special = zeta(p + 1) if b == 0 else polylog(p + 1, exp(-b))

    return prefactor*special, cond
