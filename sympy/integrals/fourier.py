"""Fourier Transforms"""
from __future__ import annotations
from sympy.core import S, pi
from sympy.core.add import Add
from sympy.core.cache import cacheit
from sympy.core.function import (
    AppliedUndef, expand, expand_mul, expand_trig)
from sympy.core.symbol import Dummy, Wild
from sympy.functions.elementary.complexes import Abs, re
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.special.delta_functions import Heaviside
from sympy.functions.special.singularity_functions import SingularityFunction
from sympy.integrals.laplace import _laplace_deep_collect
from sympy.integrals.transforms import (
    FourierTypeTransform, IntegralTransformError, _fourier_type_transform)
from sympy.logic.boolalg import And


@cacheit
def _fourier_build_rules():
    """
    This is an internal helper function that returns the table of Fourier
    transform rules in terms of the space variable `x` and the frequency
    variable `k`.  It is used by ``_fourier_apply_simple_rules``.  Each entry
    is a tuple containing:

        (space domain pattern,
         frequency domain replacement,
         condition for the rule to be applied,
         preparation function)

    The preparation function is a function with one argument that is applied
    to the expression before matching. For most rules it should be
    ``_laplace_deep_collect``.
    """
    x = Dummy('x')
    k = Dummy('k')
    a = Wild('a', exclude=[x])
    def dco(f): return _laplace_deep_collect(f, x)

    fourier_transform_rules = [
        # Gaussian
        (exp(-a*x**2), sqrt(pi)*exp(-pi**2*k**2/a)/sqrt(a),
         re(a) > 0, dco),
        # Two-sided exponential decay
        (exp(-a*Abs(x)), 2*a/(a**2 + 4*pi**2*k**2),
         re(a) > 0, dco),
        # Lorentzian
        (1/(a**2 + x**2), pi*exp(-2*pi*a*Abs(k))/a,
         re(a) > 0, dco),
    ]
    return fourier_transform_rules, x, k


def _fourier_apply_simple_rules(f, x, k):
    """
    This function applies all simple rules and returns the result if one
    of them gives a result.
    """
    simple_rules, x_, k_ = _fourier_build_rules()
    prep_old = ''
    prep_f = ''
    for x_dom, k_dom, check, prep in simple_rules:
        if prep_old != prep:
            prep_f = prep(f.subs({x: x_}))
            prep_old = prep
        ma = prep_f.match(x_dom)
        if ma:
            try:
                c = check.xreplace(ma)
            except TypeError:
                # This may happen if the space function has imaginary
                # numbers in it. Then we give up.
                continue
            if c == S.true:
                return (k_dom.xreplace(ma).subs({k_: k}), S.true)
    return None


def _fourier_expand(f, x, k):
    """
    This function tries to expand its argument with successively stronger
    methods: first it will expand on the top level, then it will expand any
    multiplications in depth, then it will try all available expansion
    methods, and finally it will try to expand trigonometric functions.

    If it can expand, it will then compute the Fourier transform of the
    expanded term.
    """

    r = expand(f, deep=False)
    if r.is_Add:
        return _fourier_transform(r, x, k, simplify=False)
    r = expand_mul(f)
    if r.is_Add:
        return _fourier_transform(r, x, k, simplify=False)
    r = expand(f)
    if r.is_Add:
        return _fourier_transform(r, x, k, simplify=False)
    if r != f:
        return _fourier_transform(r, x, k, simplify=False)
    r = expand(expand_trig(f))
    if r.is_Add:
        return _fourier_transform(r, x, k, simplify=False)
    return None


def _fourier_transform_integration(f, x, k):
    """
    The backend function for doing Fourier transforms by integration.

    This backend assumes that the frontend has already split sums such
    that `f` is not an addition anymore.  It is called with ``simplify=True``
    because this has always been the behavior of the integration backend;
    simplification of the assembled result is controlled by the
    ``simplify`` hint of the frontend instead.
    """
    try:
        return _fourier_type_transform(f, x, k, S.One, -2*S.Pi, 'Fourier',
                                       simplify=True, noconds=False)
    except IntegralTransformError:
        return None


def _fourier_transform(fn, x_, k_, *, simplify):
    """
    Front-end function of the Fourier transform. It tries to apply all known
    rules recursively, and if everything else fails, it tries to integrate.
    """
    terms_x = Add.make_args(fn)
    terms_k = []
    conditions = []

    for ff in terms_x:
        m, ft = ff.as_independent(x_, as_Add=False)
        if ft.has(SingularityFunction):
            ft = ft.rewrite(Heaviside)
        if ((r := _fourier_apply_simple_rules(ft, x_, k_)) is not None or
                (r := _fourier_expand(ft, x_, k_)) is not None):
            pass
        elif any(undef.has(x_) for undef in ft.atoms(AppliedUndef)):
            # If there are undefined functions f(x) then integration is
            # unlikely to do anything useful so we skip it and give an
            # unevaluated FourierTransform.
            r = (FourierTransform(ft, x_, k_), S.true)
        elif (r := _fourier_transform_integration(
                ft, x_, k_)) is not None:
            pass
        else:
            r = (FourierTransform(ft, x_, k_), S.true)
        (ri_, ci_) = r
        terms_k.append(m*ri_)
        conditions.append(ci_)

    result = Add(*terms_k)
    if simplify:
        result = result.simplify(doit=False)
    condition = And(*conditions)

    return result, condition


class FourierTransform(FourierTypeTransform):
    """
    Class representing unevaluated Fourier transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute Fourier transforms, see the :func:`fourier_transform`
    docstring.

    If this is called with ``.doit()``, it returns the Fourier transform as an
    expression. If it is called with ``.doit(noconds=False)``, it returns a
    tuple containing the same expression and conditions.
    """

    _name = 'Fourier'

    def a(self):
        return 1

    def b(self):
        return -2*S.Pi

    def doit(self, **hints):
        """
        Try to evaluate the transform in closed form.

        Explanation
        ===========

        Standard hints are the following:
        - ``noconds``: if True, do not return convergence conditions. The
        default setting is `True`.
        - ``simplify``: if True, it simplifies the final result. The
        default setting is `False`.
        """
        _noconds = hints.get('noconds', True)
        _simplify = hints.get('simplify', False)

        x_ = self.function_variable
        k_ = self.transform_variable
        fn = self.function

        r = _fourier_transform(fn, x_, k_, simplify=_simplify)

        if _noconds:
            return r[0]
        else:
            return r


def fourier_transform(f, x, k, **hints):
    r"""
    Compute the unitary, ordinary-frequency Fourier transform of ``f``, defined
    as

    .. math:: F(k) = \int_{-\infty}^\infty f(x) e^{-2\pi i x k} \mathrm{d} x.

    Explanation
    ===========

    If the transform cannot be computed in closed form, this
    function returns an unevaluated :class:`FourierTransform` object.

    For other Fourier transform conventions, see the function
    :func:`sympy.integrals.transforms._fourier_type_transform`.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`.
    Note that for this transform, by default ``noconds=True``.

    Examples
    ========

    >>> from sympy import fourier_transform, exp
    >>> from sympy.abc import x, k
    >>> fourier_transform(exp(-x**2), x, k)
    sqrt(pi)*exp(-pi**2*k**2)
    >>> fourier_transform(exp(-x**2), x, k, noconds=False)
    (sqrt(pi)*exp(-pi**2*k**2), True)

    See Also
    ========

    inverse_fourier_transform
    sine_transform, inverse_sine_transform
    cosine_transform, inverse_cosine_transform
    hankel_transform, inverse_hankel_transform
    mellin_transform, laplace_transform
    """
    return FourierTransform(f, x, k).doit(**hints)


class InverseFourierTransform(FourierTypeTransform):
    """
    Class representing unevaluated inverse Fourier transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute inverse Fourier transforms, see the
    :func:`inverse_fourier_transform` docstring.
    """

    _name = 'Inverse Fourier'

    def a(self):
        return 1

    def b(self):
        return 2*S.Pi


def inverse_fourier_transform(F, k, x, **hints):
    r"""
    Compute the unitary, ordinary-frequency inverse Fourier transform of `F`,
    defined as

    .. math:: f(x) = \int_{-\infty}^\infty F(k) e^{2\pi i x k} \mathrm{d} k.

    Explanation
    ===========

    If the transform cannot be computed in closed form, this
    function returns an unevaluated :class:`InverseFourierTransform` object.

    For other Fourier transform conventions, see the function
    :func:`sympy.integrals.transforms._fourier_type_transform`.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`.
    Note that for this transform, by default ``noconds=True``.

    Examples
    ========

    >>> from sympy import inverse_fourier_transform, exp, sqrt, pi
    >>> from sympy.abc import x, k
    >>> inverse_fourier_transform(sqrt(pi)*exp(-(pi*k)**2), k, x)
    exp(-x**2)
    >>> inverse_fourier_transform(sqrt(pi)*exp(-(pi*k)**2), k, x, noconds=False)
    (exp(-x**2), True)

    See Also
    ========

    fourier_transform
    sine_transform, inverse_sine_transform
    cosine_transform, inverse_cosine_transform
    hankel_transform, inverse_hankel_transform
    mellin_transform, laplace_transform
    """
    return InverseFourierTransform(F, k, x).doit(**hints)
