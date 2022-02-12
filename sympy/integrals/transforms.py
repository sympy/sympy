""" Integral Transforms """
from functools import reduce, wraps
from itertools import repeat
from sympy.core import S, pi, I
from sympy.core.add import Add
from sympy.core.function import (AppliedUndef, count_ops, Derivative, expand,
                                 expand_complex, expand_mul, Function, Lambda,
                                 WildFunction)
from sympy.core.mul import Mul
from sympy.core.numbers import igcd, ilcm
from sympy.core.relational import _canonical, Ge, Gt, Lt, Unequality, Eq
from sympy.core.sorting import default_sort_key, ordered
from sympy.core.symbol import Dummy, symbols, Wild
from sympy.core.traversal import postorder_traversal
from sympy.functions.combinatorial.factorials import factorial, rf
from sympy.functions.elementary.complexes import (re, arg, Abs, polar_lift,
                                                  periodic_argument)
from sympy.functions.elementary.exponential import exp, log, exp_polar
from sympy.functions.elementary.hyperbolic import cosh, coth, sinh, tanh, asinh
from sympy.functions.elementary.integers import ceiling
from sympy.functions.elementary.miscellaneous import Max, Min, sqrt
from sympy.functions.elementary.piecewise import Piecewise, piecewise_fold
from sympy.functions.elementary.trigonometric import cos, cot, sin, tan, atan
from sympy.functions.special.bessel import besseli, besselj, besselk, bessely
from sympy.functions.special.delta_functions import DiracDelta, Heaviside
from sympy.functions.special.error_functions import erf, erfc, Ei
from sympy.functions.special.gamma_functions import digamma, gamma, lowergamma
from sympy.functions.special.hyper import meijerg
from sympy.integrals import integrate, Integral
from sympy.integrals.meijerint import _dummy
from sympy.logic.boolalg import to_cnf, conjuncts, disjuncts, Or, And
from sympy.matrices.matrices import MatrixBase
from sympy.polys.matrices.linsolve import _lin_eq2dict, PolyNonlinearError
from sympy.polys.polyroots import roots
from sympy.polys.polytools import factor, Poly
from sympy.polys.rationaltools import together
from sympy.polys.rootoftools import CRootOf, RootSum
from sympy.simplify import simplify, hyperexpand
from sympy.simplify.powsimp import powdenest
from sympy.solvers.inequalities import _solve_inequality
from sympy.utilities.exceptions import (sympy_deprecation_warning,
                                        SymPyDeprecationWarning,
                                        ignore_warnings)
from sympy.utilities.iterables import iterable
from sympy.utilities.misc import debug


##########################################################################
# Helpers / Utilities
##########################################################################


class IntegralTransformError(NotImplementedError):
    """
    Exception raised in relation to problems computing transforms.

    Explanation
    ===========

    This class is mostly used internally; if integrals cannot be computed
    objects representing unevaluated transforms are usually returned.

    The hint ``needeval=True`` can be used to disable returning transform
    objects, and instead raise this exception if an integral cannot be
    computed.
    """
    def __init__(self, transform, function, msg):
        super().__init__(
            "%s Transform could not be computed: %s." % (transform, msg))
        self.function = function


class IntegralTransform(Function):
    """
    Base class for integral transforms.

    Explanation
    ===========

    This class represents unevaluated transforms.

    To implement a concrete transform, derive from this class and implement
    the ``_compute_transform(f, x, s, **hints)`` and ``_as_integral(f, x, s)``
    functions. If the transform cannot be computed, raise :obj:`IntegralTransformError`.

    Also set ``cls._name``. For instance,

    >>> from sympy import LaplaceTransform
    >>> LaplaceTransform._name
    'Laplace'

    Implement ``self._collapse_extra`` if your function returns more than just a
    number and possibly a convergence condition.
    """

    @property
    def function(self):
        """ The function to be transformed. """
        return self.args[0]

    @property
    def function_variable(self):
        """ The dependent variable of the function to be transformed. """
        return self.args[1]

    @property
    def transform_variable(self):
        """ The independent transform variable. """
        return self.args[2]

    @property
    def free_symbols(self):
        """
        This method returns the symbols that will exist when the transform
        is evaluated.
        """
        return self.function.free_symbols.union({self.transform_variable}) \
            - {self.function_variable}

    def _compute_transform(self, f, x, s, **hints):
        raise NotImplementedError

    def _as_integral(self, f, x, s):
        raise NotImplementedError

    def _collapse_extra(self, extra):
        cond = And(*extra)
        if cond == False:
            raise IntegralTransformError(self.__class__.name, None, '')
        return cond

    def _try_directly(self, **hints):
        T = None
        try_directly = not any(func.has(self.function_variable)
                               for func in self.function.atoms(AppliedUndef))
        if try_directly:
            try:
                T = self._compute_transform(self.function,
                    self.function_variable, self.transform_variable, **hints)
            except IntegralTransformError:
                T = None

        fn = self.function
        if not fn.is_Add:
            fn = expand_mul(fn)
        return fn, T


    def doit(self, **hints):
        """
        Try to evaluate the transform in closed form.

        Explanation
        ===========

        This general function handles linearity, but apart from that leaves
        pretty much everything to _compute_transform.

        Standard hints are the following:

        - ``simplify``: whether or not to simplify the result
        - ``noconds``: if True, do not return convergence conditions
        - ``needeval``: if True, raise IntegralTransformError instead of
                        returning IntegralTransform objects

        The default values of these hints depend on the concrete transform,
        usually the default is
        ``(simplify, noconds, needeval) = (True, False, False)``.
        """
        needeval = hints.pop('needeval', False)
        simplify = hints.pop('simplify', True)
        hints['simplify'] = simplify

        fn, T = self._try_directly(**hints)

        if T is not None:
            return T

        if fn.is_Add:
            hints['needeval'] = needeval
            res = [self.__class__(*([x] + list(self.args[1:]))).doit(**hints)
                   for x in fn.args]
            extra = []
            ress = []
            for x in res:
                if not isinstance(x, tuple):
                    x = [x]
                ress.append(x[0])
                if len(x) == 2:
                    # only a condition
                    extra.append(x[1])
                elif len(x) > 2:
                    # some region parameters and a condition (Mellin, Laplace)
                    extra += [x[1:]]
            if simplify==True:
                res = Add(*ress).simplify()
            else:
                res = Add(*ress)
            if not extra:
                return res
            try:
                extra = self._collapse_extra(extra)
                if iterable(extra):
                    return tuple([res]) + tuple(extra)
                else:
                    return (res, extra)
            except IntegralTransformError:
                pass

        if needeval:
            raise IntegralTransformError(
                self.__class__._name, self.function, 'needeval')

        # TODO handle derivatives etc

        # pull out constant coefficients
        coeff, rest = fn.as_coeff_mul(self.function_variable)
        return coeff*self.__class__(*([Mul(*rest)] + list(self.args[1:])))

    @property
    def as_integral(self):
        return self._as_integral(self.function, self.function_variable,
                                 self.transform_variable)

    def _eval_rewrite_as_Integral(self, *args, **kwargs):
        return self.as_integral


def _simplify(expr, doit):
    if doit:
        return simplify(powdenest(piecewise_fold(expr), polar=True))
    return expr


def _noconds_(default):
    """
    This is a decorator generator for dropping convergence conditions.

    Explanation
    ===========

    Suppose you define a function ``transform(*args)`` which returns a tuple of
    the form ``(result, cond1, cond2, ...)``.

    Decorating it ``@_noconds_(default)`` will add a new keyword argument
    ``noconds`` to it. If ``noconds=True``, the return value will be altered to
    be only ``result``, whereas if ``noconds=False`` the return value will not
    be altered.

    The default value of the ``noconds`` keyword will be ``default`` (i.e. the
    argument of this function).
    """
    def make_wrapper(func):
        @wraps(func)
        def wrapper(*args, noconds=default, **kwargs):
            res = func(*args, **kwargs)
            if noconds:
                return res[0]
            return res
        return wrapper
    return make_wrapper
_noconds = _noconds_(False)


##########################################################################
# Mellin Transform
##########################################################################

def _default_integrator(f, x):
    return integrate(f, (x, S.Zero, S.Infinity))


@_noconds
def _mellin_transform(f, x, s_, integrator=_default_integrator, simplify=True):
    """ Backend function to compute Mellin transforms. """
    # We use a fresh dummy, because assumptions on s might drop conditions on
    # convergence of the integral.
    s = _dummy('s', 'mellin-transform', f)
    F = integrator(x**(s - 1) * f, x)

    if not F.has(Integral):
        return _simplify(F.subs(s, s_), simplify), (S.NegativeInfinity, S.Infinity), S.true

    if not F.is_Piecewise:  # XXX can this work if integration gives continuous result now?
        raise IntegralTransformError('Mellin', f, 'could not compute integral')

    F, cond = F.args[0]
    if F.has(Integral):
        raise IntegralTransformError(
            'Mellin', f, 'integral in unexpected form')

    def process_conds(cond):
        """
        Turn ``cond`` into a strip (a, b), and auxiliary conditions.
        """
        a = S.NegativeInfinity
        b = S.Infinity
        aux = S.true
        conds = conjuncts(to_cnf(cond))
        t = Dummy('t', real=True)
        for c in conds:
            a_ = S.Infinity
            b_ = S.NegativeInfinity
            aux_ = []
            for d in disjuncts(c):
                d_ = d.replace(
                    re, lambda x: x.as_real_imag()[0]).subs(re(s), t)
                if not d.is_Relational or \
                    d.rel_op in ('==', '!=') \
                        or d_.has(s) or not d_.has(t):
                    aux_ += [d]
                    continue
                soln = _solve_inequality(d_, t)
                if not soln.is_Relational or \
                        soln.rel_op in ('==', '!='):
                    aux_ += [d]
                    continue
                if soln.lts == t:
                    b_ = Max(soln.gts, b_)
                else:
                    a_ = Min(soln.lts, a_)
            if a_ is not S.Infinity and a_ != b:
                a = Max(a_, a)
            elif b_ is not S.NegativeInfinity and b_ != a:
                b = Min(b_, b)
            else:
                aux = And(aux, Or(*aux_))
        return a, b, aux

    conds = [process_conds(c) for c in disjuncts(cond)]
    conds = [x for x in conds if x[2] != False]
    conds.sort(key=lambda x: (x[0] - x[1], count_ops(x[2])))

    if not conds:
        raise IntegralTransformError('Mellin', f, 'no convergence found')

    a, b, aux = conds[0]
    return _simplify(F.subs(s, s_), simplify), (a, b), aux


class MellinTransform(IntegralTransform):
    """
    Class representing unevaluated Mellin transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute Mellin transforms, see the :func:`mellin_transform`
    docstring.
    """

    _name = 'Mellin'

    def _compute_transform(self, f, x, s, **hints):
        return _mellin_transform(f, x, s, **hints)

    def _as_integral(self, f, x, s):
        return Integral(f*x**(s - 1), (x, S.Zero, S.Infinity))

    def _collapse_extra(self, extra):
        a = []
        b = []
        cond = []
        for (sa, sb), c in extra:
            a += [sa]
            b += [sb]
            cond += [c]
        res = (Max(*a), Min(*b)), And(*cond)
        if (res[0][0] >= res[0][1]) == True or res[1] == False:
            raise IntegralTransformError(
                'Mellin', None, 'no combined convergence.')
        return res


def mellin_transform(f, x, s, **hints):
    r"""
    Compute the Mellin transform `F(s)` of `f(x)`,

    .. math :: F(s) = \int_0^\infty x^{s-1} f(x) \mathrm{d}x.

    For all "sensible" functions, this converges absolutely in a strip
      `a < \operatorname{Re}(s) < b`.

    Explanation
    ===========

    The Mellin transform is related via change of variables to the Fourier
    transform, and also to the (bilateral) Laplace transform.

    This function returns ``(F, (a, b), cond)``
    where ``F`` is the Mellin transform of ``f``, ``(a, b)`` is the fundamental strip
    (as above), and ``cond`` are auxiliary convergence conditions.

    If the integral cannot be computed in closed form, this function returns
    an unevaluated :class:`MellinTransform` object.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`. If ``noconds=False``,
    then only `F` will be returned (i.e. not ``cond``, and also not the strip
    ``(a, b)``).

    Examples
    ========

    >>> from sympy import mellin_transform, exp
    >>> from sympy.abc import x, s
    >>> mellin_transform(exp(-x), x, s)
    (gamma(s), (0, oo), True)

    See Also
    ========

    inverse_mellin_transform, laplace_transform, fourier_transform
    hankel_transform, inverse_hankel_transform
    """
    return MellinTransform(f, x, s).doit(**hints)


def _rewrite_sin(m_n, s, a, b):
    """
    Re-write the sine function ``sin(m*s + n)`` as gamma functions, compatible
    with the strip (a, b).

    Return ``(gamma1, gamma2, fac)`` so that ``f == fac/(gamma1 * gamma2)``.

    Examples
    ========

    >>> from sympy.integrals.transforms import _rewrite_sin
    >>> from sympy import pi, S
    >>> from sympy.abc import s
    >>> _rewrite_sin((pi, 0), s, 0, 1)
    (gamma(s), gamma(1 - s), pi)
    >>> _rewrite_sin((pi, 0), s, 1, 0)
    (gamma(s - 1), gamma(2 - s), -pi)
    >>> _rewrite_sin((pi, 0), s, -1, 0)
    (gamma(s + 1), gamma(-s), -pi)
    >>> _rewrite_sin((pi, pi/2), s, S(1)/2, S(3)/2)
    (gamma(s - 1/2), gamma(3/2 - s), -pi)
    >>> _rewrite_sin((pi, pi), s, 0, 1)
    (gamma(s), gamma(1 - s), -pi)
    >>> _rewrite_sin((2*pi, 0), s, 0, S(1)/2)
    (gamma(2*s), gamma(1 - 2*s), pi)
    >>> _rewrite_sin((2*pi, 0), s, S(1)/2, 1)
    (gamma(2*s - 1), gamma(2 - 2*s), -pi)
    """
    # (This is a separate function because it is moderately complicated,
    #  and I want to doctest it.)
    # We want to use pi/sin(pi*x) = gamma(x)*gamma(1-x).
    # But there is one comlication: the gamma functions determine the
    # inegration contour in the definition of the G-function. Usually
    # it would not matter if this is slightly shifted, unless this way
    # we create an undefined function!
    # So we try to write this in such a way that the gammas are
    # eminently on the right side of the strip.
    m, n = m_n

    m = expand_mul(m/pi)
    n = expand_mul(n/pi)
    r = ceiling(-m*a - n.as_real_imag()[0])  # Don't use re(n), does not expand
    return gamma(m*s + n + r), gamma(1 - n - r - m*s), (-1)**r*pi


class MellinTransformStripError(ValueError):
    """
    Exception raised by _rewrite_gamma. Mainly for internal use.
    """
    pass


def _rewrite_gamma(f, s, a, b):
    """
    Try to rewrite the product f(s) as a product of gamma functions,
    so that the inverse Mellin transform of f can be expressed as a meijer
    G function.

    Explanation
    ===========

    Return (an, ap), (bm, bq), arg, exp, fac such that
    G((an, ap), (bm, bq), arg/z**exp)*fac is the inverse Mellin transform of f(s).

    Raises IntegralTransformError or MellinTransformStripError on failure.

    It is asserted that f has no poles in the fundamental strip designated by
    (a, b). One of a and b is allowed to be None. The fundamental strip is
    important, because it determines the inversion contour.

    This function can handle exponentials, linear factors, trigonometric
    functions.

    This is a helper function for inverse_mellin_transform that will not
    attempt any transformations on f.

    Examples
    ========

    >>> from sympy.integrals.transforms import _rewrite_gamma
    >>> from sympy.abc import s
    >>> from sympy import oo
    >>> _rewrite_gamma(s*(s+3)*(s-1), s, -oo, oo)
    (([], [-3, 0, 1]), ([-2, 1, 2], []), 1, 1, -1)
    >>> _rewrite_gamma((s-1)**2, s, -oo, oo)
    (([], [1, 1]), ([2, 2], []), 1, 1, 1)

    Importance of the fundamental strip:

    >>> _rewrite_gamma(1/s, s, 0, oo)
    (([1], []), ([], [0]), 1, 1, 1)
    >>> _rewrite_gamma(1/s, s, None, oo)
    (([1], []), ([], [0]), 1, 1, 1)
    >>> _rewrite_gamma(1/s, s, 0, None)
    (([1], []), ([], [0]), 1, 1, 1)
    >>> _rewrite_gamma(1/s, s, -oo, 0)
    (([], [1]), ([0], []), 1, 1, -1)
    >>> _rewrite_gamma(1/s, s, None, 0)
    (([], [1]), ([0], []), 1, 1, -1)
    >>> _rewrite_gamma(1/s, s, -oo, None)
    (([], [1]), ([0], []), 1, 1, -1)

    >>> _rewrite_gamma(2**(-s+3), s, -oo, oo)
    (([], []), ([], []), 1/2, 1, 8)
    """
    # Our strategy will be as follows:
    # 1) Guess a constant c such that the inversion integral should be
    #    performed wrt s'=c*s (instead of plain s). Write s for s'.
    # 2) Process all factors, rewrite them independently as gamma functions in
    #    argument s, or exponentials of s.
    # 3) Try to transform all gamma functions s.t. they have argument
    #    a+s or a-s.
    # 4) Check that the resulting G function parameters are valid.
    # 5) Combine all the exponentials.

    a_, b_ = S([a, b])

    def left(c, is_numer):
        """
        Decide whether pole at c lies to the left of the fundamental strip.
        """
        # heuristically, this is the best chance for us to solve the inequalities
        c = expand(re(c))
        if a_ is None and b_ is S.Infinity:
            return True
        if a_ is None:
            return c < b_
        if b_ is None:
            return c <= a_
        if (c >= b_) == True:
            return False
        if (c <= a_) == True:
            return True
        if is_numer:
            return None
        if a_.free_symbols or b_.free_symbols or c.free_symbols:
            return None  # XXX
            #raise IntegralTransformError('Inverse Mellin', f,
            #                     'Could not determine position of singularity %s'
            #                     ' relative to fundamental strip' % c)
        raise MellinTransformStripError('Pole inside critical strip?')

    # 1)
    s_multipliers = []
    for g in f.atoms(gamma):
        if not g.has(s):
            continue
        arg = g.args[0]
        if arg.is_Add:
            arg = arg.as_independent(s)[1]
        coeff, _ = arg.as_coeff_mul(s)
        s_multipliers += [coeff]
    for g in f.atoms(sin, cos, tan, cot):
        if not g.has(s):
            continue
        arg = g.args[0]
        if arg.is_Add:
            arg = arg.as_independent(s)[1]
        coeff, _ = arg.as_coeff_mul(s)
        s_multipliers += [coeff/pi]
    s_multipliers = [Abs(x) if x.is_extended_real else x for x in s_multipliers]
    common_coefficient = S.One
    for x in s_multipliers:
        if not x.is_Rational:
            common_coefficient = x
            break
    s_multipliers = [x/common_coefficient for x in s_multipliers]
    if not (all(x.is_Rational for x in s_multipliers) and
            common_coefficient.is_extended_real):
        raise IntegralTransformError("Gamma", None, "Nonrational multiplier")
    s_multiplier = common_coefficient/reduce(ilcm, [S(x.q)
                                             for x in s_multipliers], S.One)
    if s_multiplier == common_coefficient:
        if len(s_multipliers) == 0:
            s_multiplier = common_coefficient
        else:
            s_multiplier = common_coefficient \
                *reduce(igcd, [S(x.p) for x in s_multipliers])

    f = f.subs(s, s/s_multiplier)
    fac = S.One/s_multiplier
    exponent = S.One/s_multiplier
    if a_ is not None:
        a_ *= s_multiplier
    if b_ is not None:
        b_ *= s_multiplier

    # 2)
    numer, denom = f.as_numer_denom()
    numer = Mul.make_args(numer)
    denom = Mul.make_args(denom)
    args = list(zip(numer, repeat(True))) + list(zip(denom, repeat(False)))

    facs = []
    dfacs = []
    # *_gammas will contain pairs (a, c) representing Gamma(a*s + c)
    numer_gammas = []
    denom_gammas = []
    # exponentials will contain bases for exponentials of s
    exponentials = []

    def exception(fact):
        return IntegralTransformError("Inverse Mellin", f, "Unrecognised form '%s'." % fact)
    while args:
        fact, is_numer = args.pop()
        if is_numer:
            ugammas, lgammas = numer_gammas, denom_gammas
            ufacs = facs
        else:
            ugammas, lgammas = denom_gammas, numer_gammas
            ufacs = dfacs

        def linear_arg(arg):
            """ Test if arg is of form a*s+b, raise exception if not. """
            if not arg.is_polynomial(s):
                raise exception(fact)
            p = Poly(arg, s)
            if p.degree() != 1:
                raise exception(fact)
            return p.all_coeffs()

        # constants
        if not fact.has(s):
            ufacs += [fact]
        # exponentials
        elif fact.is_Pow or isinstance(fact, exp):
            if fact.is_Pow:
                base = fact.base
                exp_ = fact.exp
            else:
                base = exp_polar(1)
                exp_ = fact.exp
            if exp_.is_Integer:
                cond = is_numer
                if exp_ < 0:
                    cond = not cond
                args += [(base, cond)]*Abs(exp_)
                continue
            elif not base.has(s):
                a, b = linear_arg(exp_)
                if not is_numer:
                    base = 1/base
                exponentials += [base**a]
                facs += [base**b]
            else:
                raise exception(fact)
        # linear factors
        elif fact.is_polynomial(s):
            p = Poly(fact, s)
            if p.degree() != 1:
                # We completely factor the poly. For this we need the roots.
                # Now roots() only works in some cases (low degree), and CRootOf
                # only works without parameters. So try both...
                coeff = p.LT()[1]
                rs = roots(p, s)
                if len(rs) != p.degree():
                    rs = CRootOf.all_roots(p)
                ufacs += [coeff]
                args += [(s - c, is_numer) for c in rs]
                continue
            a, c = p.all_coeffs()
            ufacs += [a]
            c /= -a
            # Now need to convert s - c
            if left(c, is_numer):
                ugammas += [(S.One, -c + 1)]
                lgammas += [(S.One, -c)]
            else:
                ufacs += [-1]
                ugammas += [(S.NegativeOne, c + 1)]
                lgammas += [(S.NegativeOne, c)]
        elif isinstance(fact, gamma):
            a, b = linear_arg(fact.args[0])
            if is_numer:
                if (a > 0 and (left(-b/a, is_numer) == False)) or \
                   (a < 0 and (left(-b/a, is_numer) == True)):
                    raise NotImplementedError(
                        'Gammas partially over the strip.')
            ugammas += [(a, b)]
        elif isinstance(fact, sin):
            # We try to re-write all trigs as gammas. This is not in
            # general the best strategy, since sometimes this is impossible,
            # but rewriting as exponentials would work. However trig functions
            # in inverse mellin transforms usually all come from simplifying
            # gamma terms, so this should work.
            a = fact.args[0]
            if is_numer:
                # No problem with the poles.
                gamma1, gamma2, fac_ = gamma(a/pi), gamma(1 - a/pi), pi
            else:
                gamma1, gamma2, fac_ = _rewrite_sin(linear_arg(a), s, a_, b_)
            args += [(gamma1, not is_numer), (gamma2, not is_numer)]
            ufacs += [fac_]
        elif isinstance(fact, tan):
            a = fact.args[0]
            args += [(sin(a, evaluate=False), is_numer),
                     (sin(pi/2 - a, evaluate=False), not is_numer)]
        elif isinstance(fact, cos):
            a = fact.args[0]
            args += [(sin(pi/2 - a, evaluate=False), is_numer)]
        elif isinstance(fact, cot):
            a = fact.args[0]
            args += [(sin(pi/2 - a, evaluate=False), is_numer),
                     (sin(a, evaluate=False), not is_numer)]
        else:
            raise exception(fact)

    fac *= Mul(*facs)/Mul(*dfacs)

    # 3)
    an, ap, bm, bq = [], [], [], []
    for gammas, plus, minus, is_numer in [(numer_gammas, an, bm, True),
                                          (denom_gammas, bq, ap, False)]:
        while gammas:
            a, c = gammas.pop()
            if a != -1 and a != +1:
                # We use the gamma function multiplication theorem.
                p = Abs(S(a))
                newa = a/p
                newc = c/p
                if not a.is_Integer:
                    raise TypeError("a is not an integer")
                for k in range(p):
                    gammas += [(newa, newc + k/p)]
                if is_numer:
                    fac *= (2*pi)**((1 - p)/2) * p**(c - S.Half)
                    exponentials += [p**a]
                else:
                    fac /= (2*pi)**((1 - p)/2) * p**(c - S.Half)
                    exponentials += [p**(-a)]
                continue
            if a == +1:
                plus.append(1 - c)
            else:
                minus.append(c)

    # 4)
    # TODO

    # 5)
    arg = Mul(*exponentials)

    # for testability, sort the arguments
    an.sort(key=default_sort_key)
    ap.sort(key=default_sort_key)
    bm.sort(key=default_sort_key)
    bq.sort(key=default_sort_key)

    return (an, ap), (bm, bq), arg, exponent, fac


@_noconds_(True)
def _inverse_mellin_transform(F, s, x_, strip, as_meijerg=False):
    """ A helper for the real inverse_mellin_transform function, this one here
        assumes x to be real and positive. """
    x = _dummy('t', 'inverse-mellin-transform', F, positive=True)
    # Actually, we won't try integration at all. Instead we use the definition
    # of the Meijer G function as a fairly general inverse mellin transform.
    F = F.rewrite(gamma)
    for g in [factor(F), expand_mul(F), expand(F)]:
        if g.is_Add:
            # do all terms separately
            ress = [_inverse_mellin_transform(G, s, x, strip, as_meijerg,
                                              noconds=False)
                    for G in g.args]
            conds = [p[1] for p in ress]
            ress = [p[0] for p in ress]
            res = Add(*ress)
            if not as_meijerg:
                res = factor(res, gens=res.atoms(Heaviside))
            return res.subs(x, x_), And(*conds)

        try:
            a, b, C, e, fac = _rewrite_gamma(g, s, strip[0], strip[1])
        except IntegralTransformError:
            continue
        try:
            G = meijerg(a, b, C/x**e)
        except ValueError:
            continue
        if as_meijerg:
            h = G
        else:
            try:
                h = hyperexpand(G)
            except NotImplementedError:
                raise IntegralTransformError(
                    'Inverse Mellin', F, 'Could not calculate integral')

            if h.is_Piecewise and len(h.args) == 3:
                # XXX we break modularity here!
                h = Heaviside(x - Abs(C))*h.args[0].args[0] \
                    + Heaviside(Abs(C) - x)*h.args[1].args[0]
        # We must ensure that the integral along the line we want converges,
        # and return that value.
        # See [L], 5.2
        cond = [Abs(arg(G.argument)) < G.delta*pi]
        # Note: we allow ">=" here, this corresponds to convergence if we let
        # limits go to oo symmetrically. ">" corresponds to absolute convergence.
        cond += [And(Or(len(G.ap) != len(G.bq), 0 >= re(G.nu) + 1),
                     Abs(arg(G.argument)) == G.delta*pi)]
        cond = Or(*cond)
        if cond == False:
            raise IntegralTransformError(
                'Inverse Mellin', F, 'does not converge')
        return (h*fac).subs(x, x_), cond

    raise IntegralTransformError('Inverse Mellin', F, '')

_allowed = None


class InverseMellinTransform(IntegralTransform):
    """
    Class representing unevaluated inverse Mellin transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute inverse Mellin transforms, see the
    :func:`inverse_mellin_transform` docstring.
    """

    _name = 'Inverse Mellin'
    _none_sentinel = Dummy('None')
    _c = Dummy('c')

    def __new__(cls, F, s, x, a, b, **opts):
        if a is None:
            a = InverseMellinTransform._none_sentinel
        if b is None:
            b = InverseMellinTransform._none_sentinel
        return IntegralTransform.__new__(cls, F, s, x, a, b, **opts)

    @property
    def fundamental_strip(self):
        a, b = self.args[3], self.args[4]
        if a is InverseMellinTransform._none_sentinel:
            a = None
        if b is InverseMellinTransform._none_sentinel:
            b = None
        return a, b

    def _compute_transform(self, F, s, x, **hints):
        # IntegralTransform's doit will cause this hint to exist, but
        # InverseMellinTransform should ignore it
        hints.pop('simplify', True)
        global _allowed
        if _allowed is None:
            _allowed = {
                exp, gamma, sin, cos, tan, cot, cosh, sinh, tanh, coth,
                factorial, rf}
        for f in postorder_traversal(F):
            if f.is_Function and f.has(s) and f.func not in _allowed:
                raise IntegralTransformError('Inverse Mellin', F,
                                     'Component %s not recognised.' % f)
        strip = self.fundamental_strip
        return _inverse_mellin_transform(F, s, x, strip, **hints)

    def _as_integral(self, F, s, x):
        c = self.__class__._c
        return Integral(F*x**(-s), (s, c - S.ImaginaryUnit*S.Infinity, c +
                                    S.ImaginaryUnit*S.Infinity))/(2*S.Pi*S.ImaginaryUnit)


def inverse_mellin_transform(F, s, x, strip, **hints):
    r"""
    Compute the inverse Mellin transform of `F(s)` over the fundamental
    strip given by ``strip=(a, b)``.

    Explanation
    ===========

    This can be defined as

    .. math:: f(x) = \frac{1}{2\pi i} \int_{c - i\infty}^{c + i\infty} x^{-s} F(s) \mathrm{d}s,

    for any `c` in the fundamental strip. Under certain regularity
    conditions on `F` and/or `f`,
    this recovers `f` from its Mellin transform `F`
    (and vice versa), for positive real `x`.

    One of `a` or `b` may be passed as ``None``; a suitable `c` will be
    inferred.

    If the integral cannot be computed in closed form, this function returns
    an unevaluated :class:`InverseMellinTransform` object.

    Note that this function will assume x to be positive and real, regardless
    of the SymPy assumptions!

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`.

    Examples
    ========

    >>> from sympy import inverse_mellin_transform, oo, gamma
    >>> from sympy.abc import x, s
    >>> inverse_mellin_transform(gamma(s), s, x, (0, oo))
    exp(-x)

    The fundamental strip matters:

    >>> f = 1/(s**2 - 1)
    >>> inverse_mellin_transform(f, s, x, (-oo, -1))
    x*(1 - 1/x**2)*Heaviside(x - 1)/2
    >>> inverse_mellin_transform(f, s, x, (-1, 1))
    -x*Heaviside(1 - x)/2 - Heaviside(x - 1)/(2*x)
    >>> inverse_mellin_transform(f, s, x, (1, oo))
    (1/2 - x**2/2)*Heaviside(1 - x)/x

    See Also
    ========

    mellin_transform
    hankel_transform, inverse_hankel_transform
    """
    return InverseMellinTransform(F, s, x, strip[0], strip[1]).doit(**hints)


##########################################################################
# Laplace Transform
##########################################################################

def _simplifyconds(expr, s, a):
    r"""
    Naively simplify some conditions occurring in ``expr``, given that `\operatorname{Re}(s) > a`.

    Examples
    ========

    >>> from sympy.integrals.transforms import _simplifyconds as simp
    >>> from sympy.abc import x
    >>> from sympy import sympify as S
    >>> simp(abs(x**2) < 1, x, 1)
    False
    >>> simp(abs(x**2) < 1, x, 2)
    False
    >>> simp(abs(x**2) < 1, x, 0)
    Abs(x**2) < 1
    >>> simp(abs(1/x**2) < 1, x, 1)
    True
    >>> simp(S(1) < abs(x), x, 1)
    True
    >>> simp(S(1) < abs(1/x), x, 1)
    False

    >>> from sympy import Ne
    >>> simp(Ne(1, x**3), x, 1)
    True
    >>> simp(Ne(1, x**3), x, 2)
    True
    >>> simp(Ne(1, x**3), x, 0)
    Ne(1, x**3)
    """

    def power(ex):
        if ex == s:
            return 1
        if ex.is_Pow and ex.base == s:
            return ex.exp
        return None

    def bigger(ex1, ex2):
        """ Return True only if |ex1| > |ex2|, False only if |ex1| < |ex2|.
            Else return None. """
        if ex1.has(s) and ex2.has(s):
            return None
        if isinstance(ex1, Abs):
            ex1 = ex1.args[0]
        if isinstance(ex2, Abs):
            ex2 = ex2.args[0]
        if ex1.has(s):
            return bigger(1/ex2, 1/ex1)
        n = power(ex2)
        if n is None:
            return None
        try:
            if n > 0 and (Abs(ex1) <= Abs(a)**n) == True:
                return False
            if n < 0 and (Abs(ex1) >= Abs(a)**n) == True:
                return True
        except TypeError:
            pass

    def replie(x, y):
        """ simplify x < y """
        if not (x.is_positive or isinstance(x, Abs)) \
                or not (y.is_positive or isinstance(y, Abs)):
            return (x < y)
        r = bigger(x, y)
        if r is not None:
            return not r
        return (x < y)

    def replue(x, y):
        b = bigger(x, y)
        if b in (True, False):
            return True
        return Unequality(x, y)

    def repl(ex, *args):
        if ex in (True, False):
            return bool(ex)
        return ex.replace(*args)
    from sympy.simplify.radsimp import collect_abs
    expr = collect_abs(expr)
    expr = repl(expr, Lt, replie)
    expr = repl(expr, Gt, lambda x, y: replie(y, x))
    expr = repl(expr, Unequality, replue)
    return S(expr)

def expand_dirac_delta(expr):
    """
    Expand an expression involving DiractDelta to get it as a linear
    combination of DiracDelta functions.
    """
    return _lin_eq2dict(expr, expr.atoms(DiracDelta))


@_noconds
def _laplace_transform(f, t, s_, simplify=True):
    """ The backend function for Laplace transforms.

    This backend assumes that the frontend has already split sums
    such that `f` is to an addition anymore.
    """
    s = Dummy('s')
    a = Wild('a', exclude=[t])
    deltazero = []
    deltanonzero = []
    try:
        integratable, deltadict = expand_dirac_delta(f)
    except PolyNonlinearError:
        raise IntegralTransformError(
        'Laplace', f, 'could not expand DiracDelta expressions')
    for dirac_func, dirac_coeff in deltadict.items():
        p = dirac_func.match(DiracDelta(a*t))
        if p:
            deltazero.append(dirac_coeff.subs(t,0)/p[a])
        else:
            if dirac_func.args[0].subs(t,0).is_zero:
                raise IntegralTransformError('Laplace', f,\
                                             'not implemented yet.')
            else:
                deltanonzero.append(dirac_func*dirac_coeff)
    F = Add(integrate(exp(-s*t) * Add(integratable, *deltanonzero),
                      (t, S.Zero, S.Infinity)),
            Add(*deltazero))

    if not F.has(Integral):
        return _simplify(F.subs(s, s_), simplify), S.NegativeInfinity, S.true

    if not F.is_Piecewise:
        raise IntegralTransformError(
            'Laplace', f, 'could not compute integral')

    F, cond = F.args[0]
    if F.has(Integral):
        raise IntegralTransformError(
            'Laplace', f, 'integral in unexpected form')

    def process_conds(conds):
        """ Turn ``conds`` into a strip and auxiliary conditions. """
        a = S.NegativeInfinity
        aux = S.true
        conds = conjuncts(to_cnf(conds))
        p, q, w1, w2, w3, w4, w5 = symbols(
            'p q w1 w2 w3 w4 w5', cls=Wild, exclude=[s])
        patterns = (
            p*Abs(arg((s + w3)*q)) < w2,
            p*Abs(arg((s + w3)*q)) <= w2,
            Abs(periodic_argument((s + w3)**p*q, w1)) < w2,
            Abs(periodic_argument((s + w3)**p*q, w1)) <= w2,
            Abs(periodic_argument((polar_lift(s + w3))**p*q, w1)) < w2,
            Abs(periodic_argument((polar_lift(s + w3))**p*q, w1)) <= w2)
        for c in conds:
            a_ = S.Infinity
            aux_ = []
            for d in disjuncts(c):
                if d.is_Relational and s in d.rhs.free_symbols:
                    d = d.reversed
                if d.is_Relational and isinstance(d, (Ge, Gt)):
                    d = d.reversedsign
                for pat in patterns:
                    m = d.match(pat)
                    if m:
                        break
                if m:
                    if m[q].is_positive and m[w2]/m[p] == pi/2:
                        d = -re(s + m[w3]) < 0
                m = d.match(p - cos(w1*Abs(arg(s*w5))*w2)*Abs(s**w3)**w4 < 0)
                if not m:
                    m = d.match(
                        cos(p - Abs(periodic_argument(s**w1*w5, q))*w2)*Abs(s**w3)**w4 < 0)
                if not m:
                    m = d.match(
                        p - cos(Abs(periodic_argument(polar_lift(s)**w1*w5, q))*w2
                            )*Abs(s**w3)**w4 < 0)
                if m and all(m[wild].is_positive for wild in [w1, w2, w3, w4, w5]):
                    d = re(s) > m[p]
                d_ = d.replace(
                    re, lambda x: x.expand().as_real_imag()[0]).subs(re(s), t)
                if not d.is_Relational or \
                    d.rel_op in ('==', '!=') \
                        or d_.has(s) or not d_.has(t):
                    aux_ += [d]
                    continue
                soln = _solve_inequality(d_, t)
                if not soln.is_Relational or \
                        soln.rel_op in ('==', '!='):
                    aux_ += [d]
                    continue
                if soln.lts == t:
                    raise IntegralTransformError('Laplace', f,
                                         'convergence not in half-plane?')
                else:
                    a_ = Min(soln.lts, a_)
            if a_ is not S.Infinity:
                a = Max(a_, a)
            else:
                aux = And(aux, Or(*aux_))
        return a, aux.canonical if aux.is_Relational else aux

    conds = [process_conds(c) for c in disjuncts(cond)]
    conds2 = [x for x in conds if x[1] != False and x[0] is not S.NegativeInfinity]
    if not conds2:
        conds2 = [x for x in conds if x[1] != False]
    conds = list(ordered(conds2))

    def cnt(expr):
        if expr in (True, False):
            return 0
        return expr.count_ops()
    conds.sort(key=lambda x: (-x[0], cnt(x[1])))

    if not conds:
        raise IntegralTransformError('Laplace', f, 'no convergence found')
    a, aux = conds[0]  # XXX is [0] always the right one?

    def sbs(expr):
        return expr.subs(s, s_)
    if simplify:
        F = _simplifyconds(F, s, a)
        aux = _simplifyconds(aux, s, a)
    return _simplify(F.subs(s, s_), simplify), sbs(a), _canonical(sbs(aux))

def _laplace_deep_collect(f, t):
    """
    This is an internal helper function that traverses through the epression
    tree of `f(t)` and collects arguments. The purpose of it is that
    anything like `f(w*t-1*t-c)` will be written as `f((w-1)*t-c)` such that
    it can match `f(a*t+b)`.
    """
    func = f.func
    args = list(f.args)
    if len(f.args) == 0:
        return f
    else:
        for k in range(len(args)):
            args[k] = _laplace_deep_collect(args[k], t)
        if func.is_Add:
            return func(*args).collect(t)
        else:
            return func(*args)

def _laplace_build_rules(t, s):
    """
    This is an internal helper function that returns the table of Laplace
    transfrom rules in terms of the time variable `t` and the frequency
    variable `s`.  It is used by `_laplace_apply_rules`.
    """
    a = Wild('a', exclude=[t])
    b = Wild('b', exclude=[t])
    n = Wild('n', exclude=[t])
    tau = Wild('tau', exclude=[t])
    omega = Wild('omega', exclude=[t])
    dco = lambda f: _laplace_deep_collect(f,t)
    laplace_transform_rules = [
    # ( time domain,
    #   laplace domain,
    #   condition, convergence plane, preparation function )
    #
    # Catch constant (would otherwise be treated by 2.12)
    (a, a/s, S.true, S.Zero, dco),
    # DiracDelta rules
    (DiracDelta(a*t-b),
     exp(-s*b/a)/Abs(a),
     Or(And(a>0, b>=0), And(a<0, b<=0)), S.Zero, dco),
    (DiracDelta(a*t-b),
     S(0),
     Or(And(a<0, b>=0), And(a>0, b<=0)), S.Zero, dco),
    # Rules from http://eqworld.ipmnet.ru/en/auxiliary/inttrans/
    # 2.1
    (1,
     1/s,
     S.true, S.Zero, dco),
    # 2.2 expressed in terms of Heaviside
    (Heaviside(a*t-b),
     exp(-s*b/a)/s,
     And(a>0, b>0), S.Zero, dco),
    (Heaviside(a*t-b),
     (1-exp(-s*b/a))/s,
     And(a<0, b<0), S.Zero, dco),
    (Heaviside(a*t-b),
     1/s,
     And(a>0, b<=0), S.Zero, dco),
    (Heaviside(a*t-b),
     0,
     And(a<0, b>0), S.Zero, dco),
    # 2.3
    (t,
     1/s**2,
     S.true, S.Zero, dco),
    # 2.4
    (1/(a*t+b),
     -exp(-b/a*s)*Ei(-b/a*s)/a,
     a>0, S.Zero, dco),
    # 2.5 and 2.6 are covered by 2.11
    # 2.7
    (1/sqrt(a*t+b),
     sqrt(a*pi/s)*exp(b/a*s)*erfc(sqrt(b/a*s))/a,
     a>0, S.Zero, dco),
    # 2.8
    (sqrt(t)/(t+b),
     sqrt(pi/s)-pi*sqrt(b)*exp(b*s)*erfc(sqrt(b*s)),
     S.true, S.Zero, dco),
    # 2.9
    ((a*t+b)**(-S(3)/2),
     2*b**(-S(1)/2)-2*(pi*s/a)**(S(1)/2)*exp(b/a*s)*erfc(sqrt(b/a*s))/a,
     a>0, S.Zero, dco),
    # 2.10
    (t**(S(1)/2)*(t+a)**(-1),
     (pi/s)**(S(1)/2)-pi*a**(S(1)/2)*exp(a*s)*erfc(sqrt(a*s)),
     S.true, S.Zero, dco),
    # 2.11
    (1/(a*sqrt(t) + t**(3/2)),
     pi*a**(S(1)/2)*exp(a*s)*erfc(sqrt(a*s)),
     S.true, S.Zero, dco),
    # 2.12
    (t**n,
     gamma(n+1)/s**(n+1),
     n>-1, S.Zero, dco),
    # 2.13
    ((a*t+b)**n,
     lowergamma(n+1, b/a*s)*exp(-b/a*s)/s**(n+1)/a,
     And(n>-1, a>0), S.Zero, dco),
    # 2.14
    (t**n/(t+a),
     a**n*gamma(n+1)*lowergamma(-n,a*s),
     n>-1, S.Zero, dco),
    # 3.1
    (exp(a*t-tau),
     exp(-tau)/(s-a),
     S.true, a, dco),
    # 3.2
    (t*exp(a*t-tau),
     exp(-tau)/(s-a)**2,
     S.true, a, dco),
    # 3.3
    (t**n*exp(a*t),
     gamma(n+1)/(s-a)**(n+1),
     n>-1, a, dco),
    # 3.4 and 3.5 cannot be covered here because they are
    #     sums and only the individual sum terms will get here.
    # 3.6
    (exp(-a*t**2),
     sqrt(pi/4/a)*exp(s**2/4/a)*erfc(s/sqrt(4*a)),
     a>0, S.Zero, dco),
    # 3.7
    (t*exp(-a*t**2),
     1/(2*a)-2/sqrt(pi)/(4*a)**(S(3)/2)*s*erfc(s/sqrt(4*a)),
     S.true, S.Zero, dco),
    # 3.8
    (exp(-a/t),
     2*sqrt(a/s)*besselk(1, 2*sqrt(a*s)),
     a>=0, S.Zero, dco),
    # 3.9
    (sqrt(t)*exp(-a/t),
     S(1)/2*sqrt(pi/s**3)*(1+2*sqrt(a*s))*exp(-2*sqrt(a*s)),
     a>=0, S.Zero, dco),
    # 3.10
    (exp(-a/t)/sqrt(t),
     sqrt(pi/s)*exp(-2*sqrt(a*s)),
     a>=0, S.Zero, dco),
    # 3.11
    (exp(-a/t)/(t*sqrt(t)),
     sqrt(pi/a)*exp(-2*sqrt(a*s)),
     a>0, S.Zero, dco),
    # 3.12
    (t**n*exp(-a/t),
     2*(a/s)**((n+1)/2)*besselk(n+1, 2*sqrt(a*s)),
     a>0, S.Zero, dco),
    # 3.13
    (exp(-2*sqrt(a*t)),
     s**(-1)-sqrt(pi*a)*s**(-S(3)/2)*exp(a/s)*erfc(sqrt(a/s)),
     S.true, S.Zero, dco),
    # 3.14
    (exp(-2*sqrt(a*t))/sqrt(t),
     (pi/s)**(S(1)/2)*exp(a/s)*erfc(sqrt(a/s)),
     S.true, S.Zero, dco),
    # 4.1
    (sinh(a*t),
     a/(s**2-a**2),
     S.true, Abs(a), dco),
    # 4.2
    (sinh(a*t)**2,
     2*a**2/(s**3-4*a**2*s**2),
     S.true, Abs(2*a), dco),
    # 4.3
    (sinh(a*t)/t,
     log((s+a)/(s-a))/2,
     S.true, a, dco),
    # 4.4
    (t**n*sinh(a*t),
     gamma(n+1)/2*((s-a)**(-n-1)-(s+a)**(-n-1)),
     n>-2, Abs(a), dco),
    # 4.5
    (sinh(2*sqrt(a*t)),
     sqrt(pi*a)/s/sqrt(s)*exp(a/s),
     S.true, S.Zero, dco),
    # 4.6
    (sqrt(t)*sinh(2*sqrt(a*t)),
     pi**(S(1)/2)*s**(-S(5)/2)*(s/2+a)*exp(a/s)*erf(sqrt(a/s))-a**(S(1)/2)*s**(-2),
     S.true, S.Zero, dco),
    # 4.7
    (sinh(2*sqrt(a*t))/sqrt(t),
     pi**(S(1)/2)*s**(-S(1)/2)*exp(a/s)*erf(sqrt(a/s)),
     S.true, S.Zero, dco),
    # 4.8
    (sinh(sqrt(a*t))**2/sqrt(t),
     pi**(S(1)/2)/2*s**(-S(1)/2)*(exp(a/s)-1),
     S.true, S.Zero, dco),
    # 4.9
    (cosh(a*t),
     s/(s**2-a**2),
     S.true, Abs(a), dco),
    # 4.10
    (cosh(a*t)**2,
     (s**2-2*a**2)/(s**3-4*a**2*s**2),
     S.true, Abs(2*a), dco),
    # 4.11
    (t**n*cosh(a*t),
     gamma(n+1)/2*((s-a)**(-n-1)+(s+a)**(-n-1)),
     n>-1, Abs(a), dco),
    # 4.12
    (cosh(2*sqrt(a*t)),
     1/s+sqrt(pi*a)/s/sqrt(s)*exp(a/s)*erf(sqrt(a/s)),
     S.true, S.Zero, dco),
    # 4.13
    (sqrt(t)*cosh(2*sqrt(a*t)),
     pi**(S(1)/2)*s**(-S(5)/2)*(s/2+a)*exp(a/s),
     S.true, S.Zero, dco),
    # 4.14
    (cosh(2*sqrt(a*t))/sqrt(t),
     pi**(S(1)/2)*s**(-S(1)/2)*exp(a/s),
     S.true, S.Zero, dco),
    # 4.15
    (cosh(sqrt(a*t))**2/sqrt(t),
     pi**(S(1)/2)/2*s**(-S(1)/2)*(exp(a/s)+1),
     S.true, S.Zero, dco),
    # 5.1
    (log(a*t),
     -log(s/a+S.EulerGamma)/s,
     a>0, S.Zero, dco),
    # 5.2
    (log(1+a*t),
     -exp(s/a)/s*Ei(-s/a),
     S.true, S.Zero, dco),
    # 5.3
    (log(a*t+b),
     (log(b)-exp(s/b/a)/s*a*Ei(-s/b))/s*a,
     a>0, S.Zero, dco),
    # 5.4 is covered by 5.7
    # 5.5
    (log(t)/sqrt(t),
     -sqrt(pi/s)*(log(4*s)+S.EulerGamma),
     S.true, S.Zero, dco),
    # 5.6 is covered by 5.7
    # 5.7
    (t**n*log(t),
     gamma(n+1)*s**(-n-1)*(digamma(n+1)-log(s)),
     n>-1, S.Zero, dco),
    # 5.8
    (log(a*t)**2,
     ((log(s/a)+S.EulerGamma)**2+pi**2/6)/s,
     a>0, S.Zero, dco),
    # 5.9
    (exp(-a*t)*log(t),
     -(log(s+a)+S.EulerGamma)/(s+a),
     S.true, -a, dco),
    # 6.1
    (sin(omega*t),
     omega/(s**2+omega**2),
     S.true, S.Zero, dco),
    # 6.2
    (Abs(sin(omega*t)),
     omega/(s**2+omega**2)*coth(pi*s/2/omega),
     omega>0, S.Zero, dco),
    # 6.3 and 6.4 are covered by 1.8
    # 6.5 is covered by 1.8 together with 2.5
    # 6.6
    (sin(omega*t)/t,
     atan(omega/s),
     S.true, S.Zero, dco),
    # 6.7
    (sin(omega*t)**2/t,
     log(1+4*omega**2/s**2)/4,
     S.true, S.Zero, dco),
    # 6.8
    (sin(omega*t)**2/t**2,
     omega*atan(2*omega/s)-s*log(1+4*omega**2/s**2)/4,
     S.true, S.Zero, dco),
    # 6.9
    (sin(2*sqrt(a*t)),
      sqrt(pi*a)/s/sqrt(s)*exp(-a/s),
      a>0, S.Zero, dco),
    # 6.10
    (sin(2*sqrt(a*t))/t,
     pi*erf(sqrt(a/s)),
     a>0, S.Zero, dco),
    # 6.11
    (cos(omega*t),
     s/(s**2+omega**2),
     S.true, S.Zero, dco),
    # 6.12
    (cos(omega*t)**2,
     (s**2+2*omega**2)/(s**2+4*omega**2)/s,
     S.true, S.Zero, dco),
    # 6.13 is covered by 1.9 together with 2.5
    # 6.14 and 6.15 cannot be done with this method, the respective sum
    #       parts do not converge. Solve elsewhere if really needed.
    # 6.16
    (sqrt(t)*cos(2*sqrt(a*t)),
     sqrt(pi)/2*s**(-S(5)/2)*(s-2*a)*exp(-a/s),
     a>0, S.Zero, dco),
    # 6.17
    (cos(2*sqrt(a*t))/sqrt(t),
     sqrt(pi/s)*exp(-a/s),
     a>0, S.Zero, dco),
    # 6.18
    (sin(a*t)*sin(b*t),
     2*a*b*s/(s**2+(a+b)**2)/(s**2+(a-b)**2),
     S.true, S.Zero, dco),
    # 6.19
    (cos(a*t)*sin(b*t),
     b*(s**2-a**2+b**2)/(s**2+(a+b)**2)/(s**2+(a-b)**2),
     S.true, S.Zero, dco),
    # 6.20
    (cos(a*t)*cos(b*t),
     s*(s**2+a**2+b**2)/(s**2+(a+b)**2)/(s**2+(a-b)**2),
     S.true, S.Zero, dco),
    # 6.21
    (exp(b*t)*sin(a*t),
     a/((s-b)**2+a**2),
     S.true, b, dco),
    # 6.22
    (exp(b*t)*cos(a*t),
     (s-b)/((s-b)**2+a**2),
     S.true, b, dco),
    # 7.1
    (erf(a*t),
     exp(s**2/(2*a)**2)*erfc(s/(2*a))/s,
     a>0, S.Zero, dco),
    # 7.2
    (erf(sqrt(a*t)),
     sqrt(a)/sqrt(s+a)/s,
     a>0, S.Zero, dco),
    # 7.3
    (exp(a*t)*erf(sqrt(a*t)),
     sqrt(a)/sqrt(s)/(s-a),
     a>0, a, dco),
    # 7.4
    (erf(sqrt(a/t)/2),
     (1-exp(-sqrt(a*s)))/s,
     a>0, S.Zero, dco),
    # 7.5
    (erfc(sqrt(a*t)),
     (sqrt(s+a)-sqrt(a))/sqrt(s+a)/s,
     a>0, S.Zero, dco),
    # 7.6
    (exp(a*t)*erfc(sqrt(a*t)),
     1/(s+sqrt(a*s)),
     a>0, S.Zero, dco),
    # 7.7
    (erfc(sqrt(a/t)/2),
     exp(-sqrt(a*s))/s,
     a>0, S.Zero, dco),
    # 8.1, 8.2
    (besselj(n, a*t),
     a**n/(sqrt(s**2+a**2)*(s+sqrt(s**2+a**2))**n),
     And(a>0, n>-1), S.Zero, dco),
    # 8.3, 8.4
    (t**b*besselj(n, a*t),
     2**n/sqrt(pi)*gamma(n+S.Half)*a**n*(s**2+a**2)**(-n-S.Half),
     And(And(a>0, n>-S.Half), Eq(b, n)), S.Zero, dco),
    # 8.5
    (t**b*besselj(n, a*t),
     2**(n+1)/sqrt(pi)*gamma(n+S(3)/2)*a**n*s*(s**2+a**2)**(-n-S(3)/2),
     And(And(a>0, n>-1), Eq(b, n+1)), S.Zero, dco),
    # 8.6
    (besselj(0, 2*sqrt(a*t)),
     exp(-a/s)/s,
     a>0, S.Zero, dco),
    # 8.7, 8.8
    (t**(b)*besselj(n, 2*sqrt(a*t)),
     a**(n/2)*s**(-n-1)*exp(-a/s),
     And(And(a>0, n>-1), Eq(b, n*S.Half)), S.Zero, dco),
    # 8.9
    (besselj(0, a*sqrt(t**2+b*t)),
     exp(b*s-b*sqrt(s**2+a**2))/sqrt(s**2+a**2),
     b>0, S.Zero, dco),
    # 8.10, 8.11
    (besseli(n, a*t),
     a**n/(sqrt(s**2-a**2)*(s+sqrt(s**2-a**2))**n),
     And(a>0, n>-1), Abs(a), dco),
    # 8.12
    (t**b*besseli(n, a*t),
     2**n/sqrt(pi)*gamma(n+S.Half)*a**n*(s**2-a**2)**(-n-S.Half),
     And(And(a>0, n>-S.Half), Eq(b, n)), Abs(a), dco),
    # 8.13
    (t**b*besseli(n, a*t),
     2**(n+1)/sqrt(pi)*gamma(n+S(3)/2)*a**n*s*(s**2-a**2)**(-n-S(3)/2),
     And(And(a>0, n>-1), Eq(b, n+1)), Abs(a), dco),
    # 8.15, 8.16
    (t**(b)*besseli(n, 2*sqrt(a*t)),
     a**(n/2)*s**(-n-1)*exp(a/s),
     And(And(a>0, n>-1), Eq(b, n*S.Half)), S.Zero, dco),
    # 8.17
    (bessely(0, a*t),
     -2/pi*asinh(s/a)/sqrt(s**2+a**2),
     a>0, S.Zero, dco),
    # 8.18
    (besselk(0, a*t),
     (log(s+sqrt(s**2-a**2)))/(sqrt(s**2-a**2)),
     a>0, Abs(a), dco)
    ]
    return laplace_transform_rules

def _laplace_cr(f, a, c, **hints):
    """
    Internal helper function that will return `(f, a, c)` unless `**hints`
    contains `noconds=True`, in which case it will only return `f`.
    """
    conds = not hints.get('noconds', False)
    if conds:
        return f, a, c
    else:
        return f

def _laplace_rule_timescale(f, t, s, doit=True, **hints):
    r"""
    This internal helper function tries to apply the time-scaling rule of the
    Laplace transform and returns `None` if it cannot do it.

    Time-scaling means the following: if $F(s)$ is the Laplace transform of,
    $f(t)$, then, for any $a>0$, the Laplace transform of $f(at)$ will be
    $\frac1a F(\frac{s}{a})$. This scaling will also affect the transform's
    convergence plane.
    """
    _simplify = hints.pop('simplify', True)
    b = Wild('b', exclude=[t])
    g = WildFunction('g', nargs=1)
    k, func = f.as_independent(t, as_Add=False)
    ma1 = func.match(g)
    if ma1:
        arg = ma1[g].args[0].collect(t)
        ma2 = arg.match(b*t)
        if ma2 and ma2[b]>0:
            debug('_laplace_apply_rules match:')
            debug('      f:    %s ( %s, %s )'%(f, ma1, ma2))
            debug('      rule: amplitude and time scaling (1.1, 1.2)')
            if ma2[b]==1:
                if doit==True and not any(func.has(t) for func
                                          in ma1[g].atoms(AppliedUndef)):
                    return k*_laplace_transform(ma1[g].func(t), t, s,
                                                simplify=_simplify)
                else:
                    return k*LaplaceTransform(ma1[g].func(t), t, s, **hints)
            else:
                L = _laplace_apply_rules(ma1[g].func(t), t, s/ma2[b],
                                         doit=doit, **hints)
                try:
                    r, p, c = L
                    return (k/ma2[b]*r, p, c)
                except TypeError:
                    return k/ma2[b]*L
    return None

def _laplace_rule_heaviside(f, t, s, doit=True, **hints):
    """
    This internal helper function tries to transform a product containing the
    `Heaviside` function and returns `None` if it cannot do it.
    """
    hints.pop('simplify', True)
    a = Wild('a', exclude=[t])
    b = Wild('b', exclude=[t])
    y = Wild('y')
    g = WildFunction('g', nargs=1)
    k, func = f.as_independent(t, as_Add=False)
    ma1 = func.match(Heaviside(y)*g)
    if ma1:
        ma2 = ma1[y].match(t-a)
        ma3 = ma1[g].args[0].collect(t).match(t-b)
        if ma2 and ma2[a]>0 and ma3 and ma2[a]==ma3[b]:
            debug('_laplace_apply_rules match:')
            debug('      f:    %s ( %s, %s, %s )'%(f, ma1, ma2, ma3))
            debug('      rule: time shift (1.3)')
            L = _laplace_apply_rules(ma1[g].func(t), t, s, doit=doit, **hints)
            try:
                r, p, c = L
                return (k*exp(-ma2[a]*s)*r, p, c)
            except TypeError:
                return k*exp(-ma2[a]*s)*L
    return None


def _laplace_rule_exp(f, t, s, doit=True, **hints):
    """
    This internal helper function tries to transform a product containing the
    `exp` function and returns `None` if it cannot do it.
    """
    hints.pop('simplify', True)
    a = Wild('a', exclude=[t])

    y = Wild('y')
    z = Wild('z')
    k, func = f.as_independent(t, as_Add=False)
    ma1 = func.match(exp(y)*z)
    if ma1:
        ma2 = ma1[y].collect(t).match(a*t)
        if ma2:
            debug('_laplace_apply_rules match:')
            debug('      f:    %s ( %s, %s )'%(f, ma1, ma2))
            debug('      rule: multiply with exp (1.5)')
            L = _laplace_apply_rules(ma1[z], t, s-ma2[a], doit=doit, **hints)
            try:
                r, p, c = L
                return (r, p+ma2[a], c)
            except TypeError:
                return L
    return None

def _laplace_rule_trig(f, t, s, doit=True, **hints):
    """
    This internal helper function tries to transform a product containing a
    trigonometric function (`sin`, `cos`, `sinh`, `cosh`, ) and returns
    `None` if it cannot do it.
    """
    _simplify = hints.pop('simplify', True)
    a = Wild('a', exclude=[t])
    y = Wild('y')
    z = Wild('z')
    k, func = f.as_independent(t, as_Add=False)
    # All of the rules have a very similar form: trig(y)*z is matched, and then
    # two copies of the Laplace transform of z are shifted in the s Domain
    # and added with a weight; see rules 1.6 to 1.9 in
    # http://eqworld.ipmnet.ru/en/auxiliary/inttrans/laplace1.pdf
    # The parameters in the tuples are (fm, nu, s1, s2, sd):
    #   fm: Function to match
    #   nu: Number of the rule, for debug purposes
    #   s1: weight of the sum, 'I' for sin and '1' for all others
    #   s2: sign of the second copy of the Laplace transform of z
    #   sd: shift direction; shift along real or imaginary axis if `1` or `I`
    trigrules = [(sinh(y), '1.6',  1, -1, 1), (cosh(y), '1.7', 1, 1, 1),
                 (sin(y),  '1.8', -I, -1, I), (cos(y),  '1.9', 1, 1, I)]
    for trigrule in trigrules:
        fm, nu, s1, s2, sd = trigrule
        ma1 = func.match(fm*z)
        if ma1:
            ma2 = ma1[y].collect(t).match(a*t)
            if ma2:
                debug('_laplace_apply_rules match:')
                debug('      f:    %s ( %s, %s )'%(f, ma1, ma2))
                debug('      rule: multiply with %s (%s)'%(fm.func, nu))
                L = _laplace_apply_rules(ma1[z], t, s, doit=doit, **hints)
                try:
                    r, p, c = L
                    # The convergence plane changes only if the shift has been
                    # done along the real axis:
                    if sd==1:
                        cp_shift = Abs(ma2[a])
                    else:
                        cp_shift = 0
                    return ((s1*(r.subs(s, s-sd*ma2[a])+\
                                    s2*r.subs(s, s+sd*ma2[a]))).simplify()/2,
                            p+cp_shift, c)
                except TypeError:
                    if doit==True and _simplify==True:
                        return (s1*(L.subs(s, s-sd*ma2[a])+\
                                    s2*L.subs(s, s+sd*ma2[a]))).simplify()/2
                    else:
                        return (s1*(L.subs(s, s-sd*ma2[a])+\
                                    s2*L.subs(s, s+sd*ma2[a])))/2
    return None

def _laplace_rule_diff(f, t, s, doit=True, **hints):
    """
    This internal helper function tries to transform an expression containing
    a derivative of an undefined function and returns `None` if it cannot
    do it.
    """
    hints.pop('simplify', True)
    a = Wild('a', exclude=[t])
    y = Wild('y')
    n = Wild('n', exclude=[t])
    g = WildFunction('g', nargs=1)
    ma1 = f.match(a*Derivative(g, (t, n)))
    if ma1 and ma1[g].args[0] == t and ma1[n].is_integer:
        debug('_laplace_apply_rules match:')
        debug('      f:    %s'%(f,))
        debug('      rule: time derivative (1.11, 1.12)')
        d = []
        for k in range(ma1[n]):
            if k==0:
                y = ma1[g].func(t).subs(t, 0)
            else:
                y = Derivative(ma1[g].func(t), (t, k)).subs(t, 0)
            d.append(s**(ma1[n]-k-1)*y)
        r = s**ma1[n]*_laplace_apply_rules(ma1[g].func(t), t, s, doit=doit,
                                           **hints)
        return r - Add(*d)
    return None


def _laplace_apply_rules(f, t, s, doit=True, **hints):
    """
    Helper function for the class LaplaceTransform.

    This function does a Laplace transform based on rules and, after
    applying the rules, hands the rest over to `_laplace_transform`, which
    will attempt to integrate.

    If it is called with `doit=False`, then it will instead return
    `LaplaceTransform` objects.
    """

    k, func = f.as_independent(t, as_Add=False)

    simple_rules = _laplace_build_rules(t, s)
    for t_dom, s_dom, check, plane, prep in simple_rules:
        ma = prep(func).match(t_dom)
        if ma:
            debug('_laplace_apply_rules match:')
            debug('      f:    %s'%(func,))
            debug('      rule: %s o---o %s'%(t_dom, s_dom))
            try:
                debug('      try   %s'%(check,))
                c = check.xreplace(ma)
                debug('      check %s -> %s'%(check, c))
                if c==True:
                    return _laplace_cr(k*s_dom.xreplace(ma),
                                plane.xreplace(ma), S.true, **hints)
            except Exception:
                debug('_laplace_apply_rules did not match.')
    if f.has(DiracDelta):
        return None

    prog_rules = [_laplace_rule_timescale, _laplace_rule_heaviside,
                  _laplace_rule_exp, _laplace_rule_trig, _laplace_rule_diff]
    for p_rule in prog_rules:
        LT = p_rule(f, t, s, doit=doit, **hints)
        if LT is not None:
            return LT
    return None

class LaplaceTransform(IntegralTransform):
    """
    Class representing unevaluated Laplace transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute Laplace transforms, see the :func:`laplace_transform`
    docstring.
    """

    _name = 'Laplace'

    def _compute_transform(self, f, t, s, **hints):
        LT = _laplace_apply_rules(f, t, s, **hints)
        if LT is None:
            _simplify = hints.pop('simplify', True)
            debug('_laplace_apply_rules could not match function %s'%(f,))
            debug('    hints: %s'%(hints,))
            return _laplace_transform(f, t, s, simplify=_simplify, **hints)
        else:
            return LT

    def _as_integral(self, f, t, s):
        return Integral(f*exp(-s*t), (t, S.Zero, S.Infinity))

    def _collapse_extra(self, extra):
        conds = []
        planes = []
        for plane, cond in extra:
            conds.append(cond)
            planes.append(plane)
        cond = And(*conds)
        plane = Max(*planes)
        if cond == False:
            raise IntegralTransformError(
                'Laplace', None, 'No combined convergence.')
        return plane, cond

    def _try_directly(self, **hints):
        fn = self.function
        debug('----> _try_directly: %s'%(fn, ))
        t_ = self.function_variable
        s_ = self.transform_variable
        LT = None
        if not fn.is_Add:
            fn = expand_mul(fn)
            try:
                LT = self._compute_transform(fn, t_, s_, **hints)
            except IntegralTransformError:
                LT = None
        return fn, LT


def laplace_transform(f, t, s, legacy_matrix=True, **hints):
    r"""
    Compute the Laplace Transform `F(s)` of `f(t)`,

    .. math :: F(s) = \int_{0^{-}}^\infty e^{-st} f(t) \mathrm{d}t.

    Explanation
    ===========

    For all sensible functions, this converges absolutely in a
    half-plane

    .. math :: a < \operatorname{Re}(s)

    This function returns ``(F, a, cond)`` where ``F`` is the Laplace
    transform of ``f``, `a` is the half-plane of convergence, and `cond` are
    auxiliary convergence conditions.

    The implementation is rule-based, and if you are interested in which
    rules are applied, and whether integration is attemped, you can switch
    debug information on by setting ``sympy.SYMPY_DEBUG=True``.

    The lower bound is `0-`, meaning that this bound should be approached
    from the lower side. This is only necessary if distributions are involved.
    At present, it is only done if `f(t)` contains ``DiracDelta``, in which
    case the Laplace transform is computed implicitly as

    .. math :: F(s) = \lim_{\tau\to 0^{-}} \int_{\tau}^\infty e^{-st} f(t) \mathrm{d}t

    by applying rules.

    If the integral cannot be fully computed in closed form, this function
    returns an unevaluated :class:`LaplaceTransform` object.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`. If ``noconds=True``,
    only `F` will be returned (i.e. not ``cond``, and also not the plane ``a``).

    .. deprecated:: 1.9
        Legacy behavior for matrices where ``laplace_transform`` with
        ``noconds=False`` (the default) returns a Matrix whose elements are
        tuples. The behavior of ``laplace_transform`` for matrices will change
        in a future release of SymPy to return a tuple of the transformed
        Matrix and the convergence conditions for the matrix as a whole. Use
        ``legacy_matrix=False`` to enable the new behavior.

    Examples
    ========

    >>> from sympy import DiracDelta, exp, laplace_transform
    >>> from sympy.abc import t, s, a
    >>> laplace_transform(t**4, t, s)
    (24/s**5, 0, True)
    >>> laplace_transform(t**a, t, s)
    (gamma(a + 1)/(s*s**a), 0, re(a) > -1)
    >>> laplace_transform(DiracDelta(t)-a*exp(-a*t),t,s)
    (s/(a + s), Max(0, -a), True)

    See Also
    ========

    inverse_laplace_transform, mellin_transform, fourier_transform
    hankel_transform, inverse_hankel_transform

    """

    debug('\n***** laplace_transform(%s, %s, %s)'%(f, t, s))

    if isinstance(f, MatrixBase) and hasattr(f, 'applyfunc'):

        conds = not hints.get('noconds', False)

        if conds and legacy_matrix:
            sympy_deprecation_warning(
                """
Calling laplace_transform() on a Matrix with noconds=False (the default) is
deprecated. Either noconds=True or use legacy_matrix=False to get the new
behavior.
                """,
                deprecated_since_version="1.9",
                active_deprecations_target="deprecated-laplace-transform-matrix",
            )
            # Temporarily disable the deprecation warning for non-Expr objects
            # in Matrix
            with ignore_warnings(SymPyDeprecationWarning):
                return f.applyfunc(lambda fij: laplace_transform(fij, t, s, **hints))
        else:
            elements_trans = [laplace_transform(fij, t, s, **hints) for fij in f]
            if conds:
                elements, avals, conditions = zip(*elements_trans)
                f_laplace = type(f)(*f.shape, elements)
                return f_laplace, Max(*avals), And(*conditions)
            else:
                return type(f)(*f.shape, elements_trans)

    return LaplaceTransform(f, t, s).doit(**hints)


@_noconds_(True)
def _inverse_laplace_transform(F, s, t_, plane, simplify=True):
    """ The backend function for inverse Laplace transforms. """
    from sympy.integrals.meijerint import meijerint_inversion, _get_coeff_exp
    # There are two strategies we can try:
    # 1) Use inverse mellin transforms - related by a simple change of variables.
    # 2) Use the inversion integral.

    t = Dummy('t', real=True)

    def pw_simp(*args):
        """ Simplify a piecewise expression from hyperexpand. """
        # XXX we break modularity here!
        if len(args) != 3:
            return Piecewise(*args)
        arg = args[2].args[0].argument
        coeff, exponent = _get_coeff_exp(arg, t)
        e1 = args[0].args[0]
        e2 = args[1].args[0]
        return Heaviside(1/Abs(coeff) - t**exponent)*e1 \
            + Heaviside(t**exponent - 1/Abs(coeff))*e2

    if F.is_rational_function(s):
        F = F.apart(s)

    if F.is_Add:
        f = Add(*[_inverse_laplace_transform(X, s, t, plane, simplify)\
                     for X in F.args])
        return _simplify(f.subs(t, t_), simplify), True

    try:
        f, cond = inverse_mellin_transform(F, s, exp(-t), (None, S.Infinity),
                                           needeval=True, noconds=False)
    except IntegralTransformError:
        f = None
    if f is None:
        f = meijerint_inversion(F, s, t)
        if f is None:
            raise IntegralTransformError('Inverse Laplace', f, '')
        if f.is_Piecewise:
            f, cond = f.args[0]
            if f.has(Integral):
                raise IntegralTransformError('Inverse Laplace', f,
                                     'inversion integral of unrecognised form.')
        else:
            cond = S.true
        f = f.replace(Piecewise, pw_simp)

    if f.is_Piecewise:
        # many of the functions called below can't work with piecewise
        # (b/c it has a bool in args)
        return f.subs(t, t_), cond

    u = Dummy('u')

    def simp_heaviside(arg, H0=S.Half):
        a = arg.subs(exp(-t), u)
        if a.has(t):
            return Heaviside(arg, H0)
        rel = _solve_inequality(a > 0, u)
        if rel.lts == u:
            k = log(rel.gts)
            return Heaviside(t + k, H0)
        else:
            k = log(rel.lts)
            return Heaviside(-(t + k), H0)

    f = f.replace(Heaviside, simp_heaviside)

    def simp_exp(arg):
        return expand_complex(exp(arg))

    f = f.replace(exp, simp_exp)

    # TODO it would be nice to fix cosh and sinh ... simplify messes these
    #      exponentials up

    return _simplify(f.subs(t, t_), simplify), cond


class InverseLaplaceTransform(IntegralTransform):
    """
    Class representing unevaluated inverse Laplace transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute inverse Laplace transforms, see the
    :func:`inverse_laplace_transform` docstring.
    """

    _name = 'Inverse Laplace'
    _none_sentinel = Dummy('None')
    _c = Dummy('c')

    def __new__(cls, F, s, x, plane, **opts):
        if plane is None:
            plane = InverseLaplaceTransform._none_sentinel
        return IntegralTransform.__new__(cls, F, s, x, plane, **opts)

    @property
    def fundamental_plane(self):
        plane = self.args[3]
        if plane is InverseLaplaceTransform._none_sentinel:
            plane = None
        return plane

    def _compute_transform(self, F, s, t, **hints):
        return _inverse_laplace_transform(F, s, t, self.fundamental_plane, **hints)

    def _as_integral(self, F, s, t):
        c = self.__class__._c
        return Integral(exp(s*t)*F, (s, c - S.ImaginaryUnit*S.Infinity,
                                     c + S.ImaginaryUnit*S.Infinity))/(2*S.Pi*S.ImaginaryUnit)


def inverse_laplace_transform(F, s, t, plane=None, **hints):
    r"""
    Compute the inverse Laplace transform of `F(s)`, defined as

    .. math :: f(t) = \frac{1}{2\pi i} \int_{c-i\infty}^{c+i\infty} e^{st} F(s) \mathrm{d}s,

    for `c` so large that `F(s)` has no singularites in the
    half-plane `\operatorname{Re}(s) > c-\epsilon`.

    Explanation
    ===========

    The plane can be specified by
    argument ``plane``, but will be inferred if passed as None.

    Under certain regularity conditions, this recovers `f(t)` from its
    Laplace Transform `F(s)`, for non-negative `t`, and vice
    versa.

    If the integral cannot be computed in closed form, this function returns
    an unevaluated :class:`InverseLaplaceTransform` object.

    Note that this function will always assume `t` to be real,
    regardless of the SymPy assumption on `t`.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`.

    Examples
    ========

    >>> from sympy import inverse_laplace_transform, exp, Symbol
    >>> from sympy.abc import s, t
    >>> a = Symbol('a', positive=True)
    >>> inverse_laplace_transform(exp(-a*s)/s, s, t)
    Heaviside(-a + t)

    See Also
    ========

    laplace_transform, _fast_inverse_laplace
    hankel_transform, inverse_hankel_transform
    """
    if isinstance(F, MatrixBase) and hasattr(F, 'applyfunc'):
        return F.applyfunc(lambda Fij: inverse_laplace_transform(Fij, s, t, plane, **hints))
    return InverseLaplaceTransform(F, s, t, plane).doit(**hints)


def _fast_inverse_laplace(e, s, t):
    """Fast inverse Laplace transform of rational function including RootSum"""
    a, b, n = symbols('a, b, n', cls=Wild, exclude=[s])

    def _ilt(e):
        if not e.has(s):
            return e
        elif e.is_Add:
            return _ilt_add(e)
        elif e.is_Mul:
            return _ilt_mul(e)
        elif e.is_Pow:
            return _ilt_pow(e)
        elif isinstance(e, RootSum):
            return _ilt_rootsum(e)
        else:
            raise NotImplementedError

    def _ilt_add(e):
        return e.func(*map(_ilt, e.args))

    def _ilt_mul(e):
        coeff, expr = e.as_independent(s)
        if expr.is_Mul:
            raise NotImplementedError
        return coeff * _ilt(expr)

    def _ilt_pow(e):
        match = e.match((a*s + b)**n)
        if match is not None:
            nm, am, bm = match[n], match[a], match[b]
            if nm.is_Integer and nm < 0:
                return t**(-nm-1)*exp(-(bm/am)*t)/(am**-nm*gamma(-nm))
            if nm == 1:
                return exp(-(bm/am)*t) / am
        raise NotImplementedError

    def _ilt_rootsum(e):
        expr = e.fun.expr
        [variable] = e.fun.variables
        return RootSum(e.poly, Lambda(variable, together(_ilt(expr))))

    return _ilt(e)


##########################################################################
# Fourier Transform
##########################################################################

@_noconds_(True)
def _fourier_transform(f, x, k, a, b, name, simplify=True):
    r"""
    Compute a general Fourier-type transform

    .. math::

        F(k) = a \int_{-\infty}^{\infty} e^{bixk} f(x)\, dx.

    For suitable choice of *a* and *b*, this reduces to the standard Fourier
    and inverse Fourier transforms.
    """
    F = integrate(a*f*exp(b*S.ImaginaryUnit*x*k), (x, S.NegativeInfinity, S.Infinity))

    if not F.has(Integral):
        return _simplify(F, simplify), S.true

    integral_f = integrate(f, (x, S.NegativeInfinity, S.Infinity))
    if integral_f in (S.NegativeInfinity, S.Infinity, S.NaN) or integral_f.has(Integral):
        raise IntegralTransformError(name, f, 'function not integrable on real axis')

    if not F.is_Piecewise:
        raise IntegralTransformError(name, f, 'could not compute integral')

    F, cond = F.args[0]
    if F.has(Integral):
        raise IntegralTransformError(name, f, 'integral in unexpected form')

    return _simplify(F, simplify), cond


class FourierTypeTransform(IntegralTransform):
    """ Base class for Fourier transforms."""

    def a(self):
        raise NotImplementedError(
            "Class %s must implement a(self) but does not" % self.__class__)

    def b(self):
        raise NotImplementedError(
            "Class %s must implement b(self) but does not" % self.__class__)

    def _compute_transform(self, f, x, k, **hints):
        return _fourier_transform(f, x, k,
                                  self.a(), self.b(),
                                  self.__class__._name, **hints)

    def _as_integral(self, f, x, k):
        a = self.a()
        b = self.b()
        return Integral(a*f*exp(b*S.ImaginaryUnit*x*k), (x, S.NegativeInfinity, S.Infinity))


class FourierTransform(FourierTypeTransform):
    """
    Class representing unevaluated Fourier transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute Fourier transforms, see the :func:`fourier_transform`
    docstring.
    """

    _name = 'Fourier'

    def a(self):
        return 1

    def b(self):
        return -2*S.Pi


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
    :func:`sympy.integrals.transforms._fourier_transform`.

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
    :func:`sympy.integrals.transforms._fourier_transform`.

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


##########################################################################
# Fourier Sine and Cosine Transform
##########################################################################

@_noconds_(True)
def _sine_cosine_transform(f, x, k, a, b, K, name, simplify=True):
    """
    Compute a general sine or cosine-type transform
        F(k) = a int_0^oo b*sin(x*k) f(x) dx.
        F(k) = a int_0^oo b*cos(x*k) f(x) dx.

    For suitable choice of a and b, this reduces to the standard sine/cosine
    and inverse sine/cosine transforms.
    """
    F = integrate(a*f*K(b*x*k), (x, S.Zero, S.Infinity))

    if not F.has(Integral):
        return _simplify(F, simplify), S.true

    if not F.is_Piecewise:
        raise IntegralTransformError(name, f, 'could not compute integral')

    F, cond = F.args[0]
    if F.has(Integral):
        raise IntegralTransformError(name, f, 'integral in unexpected form')

    return _simplify(F, simplify), cond


class SineCosineTypeTransform(IntegralTransform):
    """
    Base class for sine and cosine transforms.
    Specify cls._kern.
    """

    def a(self):
        raise NotImplementedError(
            "Class %s must implement a(self) but does not" % self.__class__)

    def b(self):
        raise NotImplementedError(
            "Class %s must implement b(self) but does not" % self.__class__)


    def _compute_transform(self, f, x, k, **hints):
        return _sine_cosine_transform(f, x, k,
                                      self.a(), self.b(),
                                      self.__class__._kern,
                                      self.__class__._name, **hints)

    def _as_integral(self, f, x, k):
        a = self.a()
        b = self.b()
        K = self.__class__._kern
        return Integral(a*f*K(b*x*k), (x, S.Zero, S.Infinity))


class SineTransform(SineCosineTypeTransform):
    """
    Class representing unevaluated sine transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute sine transforms, see the :func:`sine_transform`
    docstring.
    """

    _name = 'Sine'
    _kern = sin

    def a(self):
        return sqrt(2)/sqrt(pi)

    def b(self):
        return S.One


def sine_transform(f, x, k, **hints):
    r"""
    Compute the unitary, ordinary-frequency sine transform of `f`, defined
    as

    .. math:: F(k) = \sqrt{\frac{2}{\pi}} \int_{0}^\infty f(x) \sin(2\pi x k) \mathrm{d} x.

    Explanation
    ===========

    If the transform cannot be computed in closed form, this
    function returns an unevaluated :class:`SineTransform` object.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`.
    Note that for this transform, by default ``noconds=True``.

    Examples
    ========

    >>> from sympy import sine_transform, exp
    >>> from sympy.abc import x, k, a
    >>> sine_transform(x*exp(-a*x**2), x, k)
    sqrt(2)*k*exp(-k**2/(4*a))/(4*a**(3/2))
    >>> sine_transform(x**(-a), x, k)
    2**(1/2 - a)*k**(a - 1)*gamma(1 - a/2)/gamma(a/2 + 1/2)

    See Also
    ========

    fourier_transform, inverse_fourier_transform
    inverse_sine_transform
    cosine_transform, inverse_cosine_transform
    hankel_transform, inverse_hankel_transform
    mellin_transform, laplace_transform
    """
    return SineTransform(f, x, k).doit(**hints)


class InverseSineTransform(SineCosineTypeTransform):
    """
    Class representing unevaluated inverse sine transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute inverse sine transforms, see the
    :func:`inverse_sine_transform` docstring.
    """

    _name = 'Inverse Sine'
    _kern = sin

    def a(self):
        return sqrt(2)/sqrt(pi)

    def b(self):
        return S.One


def inverse_sine_transform(F, k, x, **hints):
    r"""
    Compute the unitary, ordinary-frequency inverse sine transform of `F`,
    defined as

    .. math:: f(x) = \sqrt{\frac{2}{\pi}} \int_{0}^\infty F(k) \sin(2\pi x k) \mathrm{d} k.

    Explanation
    ===========

    If the transform cannot be computed in closed form, this
    function returns an unevaluated :class:`InverseSineTransform` object.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`.
    Note that for this transform, by default ``noconds=True``.

    Examples
    ========

    >>> from sympy import inverse_sine_transform, exp, sqrt, gamma
    >>> from sympy.abc import x, k, a
    >>> inverse_sine_transform(2**((1-2*a)/2)*k**(a - 1)*
    ...     gamma(-a/2 + 1)/gamma((a+1)/2), k, x)
    x**(-a)
    >>> inverse_sine_transform(sqrt(2)*k*exp(-k**2/(4*a))/(4*sqrt(a)**3), k, x)
    x*exp(-a*x**2)

    See Also
    ========

    fourier_transform, inverse_fourier_transform
    sine_transform
    cosine_transform, inverse_cosine_transform
    hankel_transform, inverse_hankel_transform
    mellin_transform, laplace_transform
    """
    return InverseSineTransform(F, k, x).doit(**hints)


class CosineTransform(SineCosineTypeTransform):
    """
    Class representing unevaluated cosine transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute cosine transforms, see the :func:`cosine_transform`
    docstring.
    """

    _name = 'Cosine'
    _kern = cos

    def a(self):
        return sqrt(2)/sqrt(pi)

    def b(self):
        return S.One


def cosine_transform(f, x, k, **hints):
    r"""
    Compute the unitary, ordinary-frequency cosine transform of `f`, defined
    as

    .. math:: F(k) = \sqrt{\frac{2}{\pi}} \int_{0}^\infty f(x) \cos(2\pi x k) \mathrm{d} x.

    Explanation
    ===========

    If the transform cannot be computed in closed form, this
    function returns an unevaluated :class:`CosineTransform` object.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`.
    Note that for this transform, by default ``noconds=True``.

    Examples
    ========

    >>> from sympy import cosine_transform, exp, sqrt, cos
    >>> from sympy.abc import x, k, a
    >>> cosine_transform(exp(-a*x), x, k)
    sqrt(2)*a/(sqrt(pi)*(a**2 + k**2))
    >>> cosine_transform(exp(-a*sqrt(x))*cos(a*sqrt(x)), x, k)
    a*exp(-a**2/(2*k))/(2*k**(3/2))

    See Also
    ========

    fourier_transform, inverse_fourier_transform,
    sine_transform, inverse_sine_transform
    inverse_cosine_transform
    hankel_transform, inverse_hankel_transform
    mellin_transform, laplace_transform
    """
    return CosineTransform(f, x, k).doit(**hints)


class InverseCosineTransform(SineCosineTypeTransform):
    """
    Class representing unevaluated inverse cosine transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute inverse cosine transforms, see the
    :func:`inverse_cosine_transform` docstring.
    """

    _name = 'Inverse Cosine'
    _kern = cos

    def a(self):
        return sqrt(2)/sqrt(pi)

    def b(self):
        return S.One


def inverse_cosine_transform(F, k, x, **hints):
    r"""
    Compute the unitary, ordinary-frequency inverse cosine transform of `F`,
    defined as

    .. math:: f(x) = \sqrt{\frac{2}{\pi}} \int_{0}^\infty F(k) \cos(2\pi x k) \mathrm{d} k.

    Explanation
    ===========

    If the transform cannot be computed in closed form, this
    function returns an unevaluated :class:`InverseCosineTransform` object.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`.
    Note that for this transform, by default ``noconds=True``.

    Examples
    ========

    >>> from sympy import inverse_cosine_transform, sqrt, pi
    >>> from sympy.abc import x, k, a
    >>> inverse_cosine_transform(sqrt(2)*a/(sqrt(pi)*(a**2 + k**2)), k, x)
    exp(-a*x)
    >>> inverse_cosine_transform(1/sqrt(k), k, x)
    1/sqrt(x)

    See Also
    ========

    fourier_transform, inverse_fourier_transform,
    sine_transform, inverse_sine_transform
    cosine_transform
    hankel_transform, inverse_hankel_transform
    mellin_transform, laplace_transform
    """
    return InverseCosineTransform(F, k, x).doit(**hints)


##########################################################################
# Hankel Transform
##########################################################################

@_noconds_(True)
def _hankel_transform(f, r, k, nu, name, simplify=True):
    r"""
    Compute a general Hankel transform

    .. math:: F_\nu(k) = \int_{0}^\infty f(r) J_\nu(k r) r \mathrm{d} r.
    """
    F = integrate(f*besselj(nu, k*r)*r, (r, S.Zero, S.Infinity))

    if not F.has(Integral):
        return _simplify(F, simplify), S.true

    if not F.is_Piecewise:
        raise IntegralTransformError(name, f, 'could not compute integral')

    F, cond = F.args[0]
    if F.has(Integral):
        raise IntegralTransformError(name, f, 'integral in unexpected form')

    return _simplify(F, simplify), cond


class HankelTypeTransform(IntegralTransform):
    """
    Base class for Hankel transforms.
    """

    def doit(self, **hints):
        return self._compute_transform(self.function,
                                       self.function_variable,
                                       self.transform_variable,
                                       self.args[3],
                                       **hints)

    def _compute_transform(self, f, r, k, nu, **hints):
        return _hankel_transform(f, r, k, nu, self._name, **hints)

    def _as_integral(self, f, r, k, nu):
        return Integral(f*besselj(nu, k*r)*r, (r, S.Zero, S.Infinity))

    @property
    def as_integral(self):
        return self._as_integral(self.function,
                                 self.function_variable,
                                 self.transform_variable,
                                 self.args[3])


class HankelTransform(HankelTypeTransform):
    """
    Class representing unevaluated Hankel transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute Hankel transforms, see the :func:`hankel_transform`
    docstring.
    """

    _name = 'Hankel'


def hankel_transform(f, r, k, nu, **hints):
    r"""
    Compute the Hankel transform of `f`, defined as

    .. math:: F_\nu(k) = \int_{0}^\infty f(r) J_\nu(k r) r \mathrm{d} r.

    Explanation
    ===========

    If the transform cannot be computed in closed form, this
    function returns an unevaluated :class:`HankelTransform` object.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`.
    Note that for this transform, by default ``noconds=True``.

    Examples
    ========

    >>> from sympy import hankel_transform, inverse_hankel_transform
    >>> from sympy import exp
    >>> from sympy.abc import r, k, m, nu, a

    >>> ht = hankel_transform(1/r**m, r, k, nu)
    >>> ht
    2*k**(m - 2)*gamma(-m/2 + nu/2 + 1)/(2**m*gamma(m/2 + nu/2))

    >>> inverse_hankel_transform(ht, k, r, nu)
    r**(-m)

    >>> ht = hankel_transform(exp(-a*r), r, k, 0)
    >>> ht
    a/(k**3*(a**2/k**2 + 1)**(3/2))

    >>> inverse_hankel_transform(ht, k, r, 0)
    exp(-a*r)

    See Also
    ========

    fourier_transform, inverse_fourier_transform
    sine_transform, inverse_sine_transform
    cosine_transform, inverse_cosine_transform
    inverse_hankel_transform
    mellin_transform, laplace_transform
    """
    return HankelTransform(f, r, k, nu).doit(**hints)


class InverseHankelTransform(HankelTypeTransform):
    """
    Class representing unevaluated inverse Hankel transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute inverse Hankel transforms, see the
    :func:`inverse_hankel_transform` docstring.
    """

    _name = 'Inverse Hankel'


def inverse_hankel_transform(F, k, r, nu, **hints):
    r"""
    Compute the inverse Hankel transform of `F` defined as

    .. math:: f(r) = \int_{0}^\infty F_\nu(k) J_\nu(k r) k \mathrm{d} k.

    Explanation
    ===========

    If the transform cannot be computed in closed form, this
    function returns an unevaluated :class:`InverseHankelTransform` object.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`.
    Note that for this transform, by default ``noconds=True``.

    Examples
    ========

    >>> from sympy import hankel_transform, inverse_hankel_transform
    >>> from sympy import exp
    >>> from sympy.abc import r, k, m, nu, a

    >>> ht = hankel_transform(1/r**m, r, k, nu)
    >>> ht
    2*k**(m - 2)*gamma(-m/2 + nu/2 + 1)/(2**m*gamma(m/2 + nu/2))

    >>> inverse_hankel_transform(ht, k, r, nu)
    r**(-m)

    >>> ht = hankel_transform(exp(-a*r), r, k, 0)
    >>> ht
    a/(k**3*(a**2/k**2 + 1)**(3/2))

    >>> inverse_hankel_transform(ht, k, r, 0)
    exp(-a*r)

    See Also
    ========

    fourier_transform, inverse_fourier_transform
    sine_transform, inverse_sine_transform
    cosine_transform, inverse_cosine_transform
    hankel_transform
    mellin_transform, laplace_transform
    """
    return InverseHankelTransform(F, k, r, nu).doit(**hints)
