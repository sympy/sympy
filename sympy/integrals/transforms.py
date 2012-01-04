""" Integral Transforms """
from sympy.integrals import integrate, Integral
from sympy.core.numbers import oo
from sympy.core.symbol import Dummy
from sympy.core.function import Function
from sympy.logic.boolalg import to_cnf, conjuncts, disjuncts, Or, And
from sympy.simplify import simplify
from sympy.core import S

from sympy.integrals.meijerint import _dummy

##########################################################################
# Helpers / Utilities
##########################################################################

class IntegralTransformError(NotImplementedError):
    """
    Exception raised in relation to problems computing transforms.

    This class is mostly used internally; if integrals cannot be computed
    objects representing unevaluated transforms are usually returned.

    The hint ``needeval=True`` can be used to disable returning transform
    objects, and instead raise this exception if an integral cannot be
    computed.
    """
    def __init__(self, transform, function, msg):
        super(IntegralTransformError, self).__init__(
            "%s Transform could not be computed: %s." % (transform, msg))
        self.function = function

class IntegralTransform(Function):
    """
    Base class for integral transforms.

    This class represents unevaluated transforms.

    To implement a concrete transform, derive from this class and implement
    the _compute_transform(f, x, s, **hints) and _as_integral(f, x, s)
    functions. If the transform cannot be computed, raise IntegralTransformError.

    Also set cls._name.

    Implement self._collapse_extra if your function returns more than just a
    number and possibly a convergence condition.
    """

    nargs = 3

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
        return self.function.free_symbols.union(set([self.transform_variable])) \
               - set([self.function_variable])

    def _compute_transform(self, f, x, s, **hints):
        raise NotImplementedError

    def _as_integral(self, f, x, s):
        raise NotImplementedError

    def _collapse_extra(self, extra):
        from sympy import And
        cond = And(*extra)
        if cond is False:
            raise IntegralTransformError(self.__class__.name, None, '')

    def doit(self, **hints):
        """
        Try to evaluate the transform in closed form.

        This general function handles linearity, but apart from that leaves
        pretty much everything to _compute_transform.

        Standard hints are the following:

        - ``simplify``: whether or not to simplify the result
        - ``noconds``: if True, don't return convergence conditions
        - ``needeval``: if True, raise IntegralTransformError instead of
                        returning IntegralTransform objects

        The default values of these hints depend on the concrete transform,
        usually the default is
        ``(simplify, noconds, needeval) = (True, False, False)``.
        """
        from sympy import Add
        from sympy.core.function import AppliedUndef
        needeval = hints.pop('needeval', False)
        try_directly = not any(func.has(self.function_variable) \
                               for func in self.function.atoms(AppliedUndef))
        if try_directly:
            try:
                return self._compute_transform(self.function,
                    self.function_variable, self.transform_variable, **hints)
            except IntegralTransformError:
                pass

        if self.function.is_Add:
            hints['needeval'] = needeval
            res = [self.__class__(*([x] + list(self.args[1:]))).doit(**hints)
                   for x in self.function.args]
            extra = []
            ress = []
            for x in res:
                if not isinstance(x, tuple):
                    x = [x]
                ress.append(x[0])
                if len(x) > 1:
                    extra += [x[1:]]
            res = Add(*ress)
            if not extra:
                return res
            try:
                extra = self._collapse_extra(extra)
                return tuple([res]) + tuple(extra)
            except IntegralTransformError:
                pass

        if needeval:
            raise IntegralTransformError(self.__class__._name, self.function, 'needeval')

        # TODO handle derivatives etc

        return self

    @property
    def as_integral(self):
        return self._as_integral(self.function, self.function_variable,
                                 self.transform_variable)

    def _eval_rewrite_as_Integral(self, *args):
        return self.as_integral

from sympy.solvers.inequalities import _solve_inequality

def _simplify(expr, doit):
    from sympy import powdenest, powsimp
    if doit:
        return simplify(powdenest(expr, polar=True))
    return expr

def _noconds_(default):
    """
    This is a decorator generator for dropping convergence conditions.

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
        from sympy.core.decorators import wraps
        @wraps(func)
        def wrapper(*args, **kwargs):
            noconds = kwargs.pop('noconds', default)
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
    return integrate(f, (x, 0, oo))

@_noconds
def _mellin_transform(f, x, s_, integrator=_default_integrator, simplify=True):
    """ Backend function to compute mellin transforms. """
    from sympy import re, Max, Min
    # We use a fresh dummy, because assumptions on s might drop conditions on
    # convergence of the integral.
    s = _dummy('s', 'mellin-transform', f)
    F = integrator(x**(s-1) * f, x)

    if not F.has(Integral):
        return _simplify(F.subs(s, s_), simplify), (-oo, oo), True

    if not F.is_Piecewise:
        raise IntegralTransformError('Mellin', f, 'could not compute integral')

    F, cond = F.args[0]
    if F.has(Integral):
        raise IntegralTransformError('Mellin', f, 'integral in unexpected form')

    a = -oo
    b = oo
    aux = True
    conds = conjuncts(to_cnf(cond))
    t = Dummy('t', real=True)
    for c in conds:
        a_ = oo
        b_ = -oo
        aux_ = []
        for d in disjuncts(c):
            d_ = d.replace(re, lambda x: x.as_real_imag()[0]).subs(re(s), t)
            if not d.is_Relational or (d.rel_op != '<' and d.rel_op != '<=') \
               or d_.has(s) or not d_.has(t):
                aux_ += [d]
                continue
            soln = _solve_inequality(d_, t)
            if not soln.is_Relational or \
               (soln.rel_op != '<' and soln.rel_op != '<='):
                aux_ += [d]
                continue
            if soln.lhs == t:
                b_ = Max(soln.rhs, b_)
            else:
                a_ = Min(soln.lhs, a_)
        if a_ != oo and a_ != b:
            a = Max(a_, a)
        elif b_ != -oo and b_ != a:
            b = Min(b_, b)
        else:
            aux = And(aux, Or(*aux_))

    if aux is False:
        raise IntegralTransformError('Mellin', f, 'no convergence found')

    return _simplify(F.subs(s, s_), simplify), (a, b), aux

class MellinTransform(IntegralTransform):
    """
    Class representing unevaluated mellin transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute mellin transforms, see the :func:`mellin_transform`
    docstring.
    """

    _name = 'Mellin'

    def _compute_transform(self, f, x, s, **hints):
        return _mellin_transform(f, x, s, **hints)

    def _as_integral(self, f, x, s):
        from sympy import Integral
        return Integral(f*x**(s-1), (x, 0, oo))

    def _collapse_extra(self, extra):
        from sympy import And, Max, Min
        a = []
        b = []
        cond = []
        for (sa, sb), c in extra:
            a += [sa]
            b += [sb]
            cond += [c]
        res = (Max(*a), Min(*b)), And(*cond)
        if (res[0][0] >= res[0][1]) is True or res[1] is False:
            raise IntegralTransformError('Mellin', None, 'no combined convergence.')
        return res

def mellin_transform(f, x, s, **hints):
    r"""
    Compute the mellin transform `F(s)` of `f(x)`,

    .. math :: F(s) = \int_0^\infty x^{s-1} f(x) \mathrm{d}x.

    For all "sensible" functions, this converges absolutely in a strip
      `a < Re(s) < b`.

    The mellin transform is related via change of variables to the fourier
    transform, and also to the (bilateral) laplace transform.

    This function returns (F, (a, b), cond)
    where `F` is the mellin transform of `f`, `(a, b)` is the fundamental strip
    (as above), and cond are auxiliary convergence conditions.

    If the integral cannot be computed in closed form, this function returns
    an unevaluated MellinTransform object.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.transforms.IntegralTransform.doit`. If ``noconds=False``,
    then only `F` will be returned (i.e. not ``cond``, and also not the strip
    ``(a, b)``).

    >>> from sympy.integrals.transforms import mellin_transform
    >>> from sympy import exp
    >>> from sympy.abc import x, s
    >>> mellin_transform(exp(-x), x, s)
    (gamma(s), (0, oo), True)

    See Also
    ========

    inverse_mellin_transform, laplace_transform, fourier_transform
    """
    return MellinTransform(f, x, s).doit(**hints)

def _rewrite_sin((m, n), s, a, b):
    """
    Re-write the sine function sin(m*s + n) as gamma functions, compatible
    with the strip (a, b).

    Return (gamma1, gamma2, fac) so that f == fac/(gamma1 * gamma2).

    >>> from sympy.integrals.transforms import _rewrite_sin
    >>> from sympy import pi, S
    >>> from sympy.abc import s
    >>> _rewrite_sin((pi, 0), s, 0, 1)
    (gamma(s), gamma(-s + 1), pi)
    >>> _rewrite_sin((pi, 0), s, 1, 0)
    (gamma(s - 1), gamma(-s + 2), -pi)
    >>> _rewrite_sin((pi, 0), s, -1, 0)
    (gamma(s + 1), gamma(-s), -pi)
    >>> _rewrite_sin((pi, pi/2), s, S(1)/2, S(3)/2)
    (gamma(s - 1/2), gamma(-s + 3/2), -pi)
    >>> _rewrite_sin((pi, pi), s, 0, 1)
    (gamma(s), gamma(-s + 1), -pi)
    >>> _rewrite_sin((2*pi, 0), s, 0, S(1)/2)
    (gamma(2*s), gamma(-2*s + 1), pi)
    >>> _rewrite_sin((2*pi, 0), s, S(1)/2, 1)
    (gamma(2*s - 1), gamma(-2*s + 2), -pi)
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
    from sympy import expand_mul, pi, ceiling, gamma, re
    m = expand_mul(m/pi)
    n = expand_mul(n/pi)
    r = ceiling(-m*a - n.as_real_imag()[0]) # Don't use re(n), does not expand
    return gamma(m*s + n + r), gamma(1 - n - r - m*s), (-1)**r*pi
def _rewrite_gamma(f, s, a, b):
    """
    Try to rewrite the product f(s) as a product of gamma functions,
    so that the inverse mellin transform of f can be expressed as a meijer
    G function.

    Return (an, ap), (bm, bq), arg, exp, fac such that
    G((an, ap), (bm, bq), arg/z**exp)*fac is the inverse mellin transform of f(s).

    Raises IntegralTransformError or ValueError on failure.

    It is asserted that f has no poles in the fundamental strip designated by
    (a, b). One of a and b is allowed to be None. The fundamental strip is
    important, because it determines the inversion contour.

    This function can handle exponentials, linear factors, trigonometric
    functions.

    This is a helper function for inverse_mellin_transform that will not
    attempt any transformations on f.

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
    from itertools import repeat
    from sympy import (Poly, gamma, Mul, re, RootOf, exp as exp_, E, expand,
                       roots, ilcm, pi, sin, cos, tan, cot, igcd)
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
        if a_ is None:
            return c < b_
        if b_ is None:
            return c <= a_
        if (c >= b_) is True:
            return False
        if (c <= a_) is True:
            return True
        if is_numer:
            return None
        if a_.free_symbols or b_.free_symbols or c.free_symbols:
            return None # XXX
            #raise IntegralTransformError('Inverse Mellin', f,
            #                     'Could not determine position of singularity %s'
            #                     ' relative to fundamental strip' % c)
        raise ValueError('Pole inside critical strip?')

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
    s_multipliers = [abs(x) for x in s_multipliers if x.is_real]
    common_coefficient = S(1)
    for x in s_multipliers:
        if not x.is_Rational:
            common_coefficient = x
            break
    s_multipliers = [x/common_coefficient for x in s_multipliers]
    if any(not x.is_Rational for x in s_multipliers):
        raise NotImplementedError
    s_multiplier = common_coefficient/reduce(ilcm, [S(x.q) for x in s_multipliers], S(1))
    if s_multiplier == common_coefficient:
        if len(s_multipliers) == 0:
            s_multiplier = common_coefficient
        else:
            s_multiplier = common_coefficient \
                           *reduce(igcd, [S(x.p) for x in s_multipliers])

    exponent = S(1)
    fac = S(1)
    f = f.subs(s, s/s_multiplier)
    fac /= s_multiplier
    exponent = 1/s_multiplier
    if a_ is not None:
        a_ *= s_multiplier
    if b_ is not None:
        b_ *= s_multiplier

    # 2)
    numer, denom = f.as_numer_denom()
    numer = Mul.make_args(numer)
    denom = Mul.make_args(denom)
    args = zip(numer, repeat(True)) + zip(denom, repeat(False))

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
            ufacs, lfacs = facs, dfacs
        else:
            ugammas, lgammas = denom_gammas, numer_gammas
            ufacs, lfacs = dfacs, facs

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
        elif fact.is_Pow or isinstance(fact, exp_):
            if fact.is_Pow:
                base = fact.base
                exp  = fact.exp
            else:
                base = E
                exp  = fact.args[0]
            if exp.is_Integer:
                cond = is_numer
                if exp < 0:
                    cond = not cond
                args += [(base, cond)]*abs(exp)
                continue
            elif not base.has(s):
                a, b = linear_arg(exp)
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
                # Now roots() only works in some cases (low degree), and RootOf
                # only works without parameters. So try both...
                coeff = p.LT()[1]
                rs = roots(p, s)
                if len(rs) != p.degree():
                    rs = RootOf.all_roots(p)
                ufacs += [coeff]
                args += [(s - c, is_numer) for c in rs]
                continue
            a, c = p.all_coeffs()
            ufacs += [a]
            c /= -a
            # Now need to convert s - c
            if left(c, is_numer):
                ugammas += [(S(1), -c + 1)]
                lgammas += [(S(1), -c)]
            else:
                ufacs += [-1]
                ugammas += [(S(-1), c + 1)]
                lgammas += [(S(-1), c)]
        elif isinstance(fact, gamma):
            a, b = linear_arg(fact.args[0])
            if is_numer:
                if (a > 0 and (left(-b/a, is_numer) is False)) or \
                   (a < 0 and (left(-b/a, is_numer) is True)):
                    raise NotImplementedError('Gammas partially over the strip.')
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
                p = abs(S(a))
                newa = a/p
                newc = c/p
                assert a.is_Integer
                for k in range(p):
                    gammas += [(newa, newc + k/p)]
                if is_numer:
                    fac *= (2*pi)**((1 - p)/2) * p**(c - S(1)/2)
                    exponentials += [p**a]
                else:
                    fac /= (2*pi)**((1 - p)/2) * p**(c - S(1)/2)
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
    an.sort()
    ap.sort()
    bm.sort()
    bq.sort()

    return (an, ap), (bm, bq), arg, exponent, fac

@_noconds_(True)
def _inverse_mellin_transform(F, s, x_, strip, as_meijerg=False):
    """ A helper for the real inverse_mellin_transform function, this one here
        assumes x to be real and positive. """
    from sympy import (expand, expand_mul, hyperexpand, meijerg, And, Or,
                       arg, pi, re, factor, Heaviside, gamma, Add)
    x = _dummy('t', 'inverse-mellin-transform', F, positive=True)
    # Actually, we won't try integration at all. Instead we use the definition
    # of the Meijer G function as a fairly general inverse mellin transform.
    F = F.rewrite(gamma)
    for g in [factor(F), expand_mul(F), expand(F)]:
        if g.is_Add:
            # do all terms separately
            ress = [_inverse_mellin_transform(G, s, x, strip, as_meijerg,
                                              noconds=False) \
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
        G = meijerg(a, b, C/x**e)
        if as_meijerg:
            h = G
        else:
            h = hyperexpand(G)
            if h.is_Piecewise and len(h.args) == 3:
                # XXX we break modularity here!
                h = Heaviside(x - abs(C))*h.args[0].args[0] \
                  + Heaviside(abs(C) - x)*h.args[1].args[0]
        # We must ensure that the intgral along the line we want converges,
        # and return that value.
        # See [L], 5.2
        cond = [abs(arg(G.argument)) < G.delta*pi]
        # Note: we allow ">=" here, this corresponds to convergence if we let
        # limits go to oo symetrically. ">" corresponds to absolute convergence.
        cond += [And(Or(len(G.ap) != len(G.bq), 0 >= re(G.nu) + 1),
                     abs(arg(G.argument)) == G.delta*pi)]
        cond = Or(*cond)
        if cond is False:
            raise IntegralTransformError('Inverse Mellin', F, 'does not converge')
        return (h*fac).subs(x, x_), cond

    raise IntegralTransformError('Inverse Mellin', F, '')

_allowed = None
class InverseMellinTransform(IntegralTransform):
    """
    Class representing unevaluated inverse mellin transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute inverse mellin transforms, see the
    :func:`inverse_mellin_transform` docstring.
    """

    nargs = 5

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
        a, b  = self.args[3], self.args[4]
        if a is InverseMellinTransform._none_sentinel:
            a = None
        if b is InverseMellinTransform._none_sentinel:
            b = None
        return a, b

    def _compute_transform(self, F, s, x, **hints):
        from sympy import postorder_traversal
        global _allowed
        if _allowed is None:
            from sympy import (exp, gamma, sin, cos, tan, cot, cosh, sinh, tanh,
                               coth, factorial, rf)
            _allowed = set([exp, gamma, sin, cos, tan, cot, cosh, sinh, tanh, coth,
                            factorial, rf])
        for f in postorder_traversal(F):
            if f.is_Function and f.has(s) and f.func not in _allowed:
                raise IntegralTransformError('Inverse Mellin', F,
                                     'Component %s not recognised.' % f)
        strip = self.fundamental_strip
        return _inverse_mellin_transform(F, s, x, strip, **hints)

    def _as_integral(self, F, s, x):
        from sympy import Integral, I, oo
        c = self.__class__._c
        return Integral(F*x**(-s), (s, c - I*oo, c + I*oo))

def inverse_mellin_transform(F, s, x, strip, **hints):
    r"""
    Compute the inverse mellin transform of `F(s)` over the fundamental
    strip given by ``strip=(a, b)``.

    This can be defined as

    .. math:: f(x) = \int_{c - i\infty}^{c + i\infty} x^{-s} F(s) \mathrm{d}s,

    for any `c` in the fundamental strip. Under certain regularity
    conditions on `F` and/or `f`,
    this recovers `f` from its mellin transform `F`
    (and vice versa), for positive real `x`.

    One of `a` or `b` may be passed as None; a suitable `c` will be
    inferred.

    If the integral cannot be computed in closed form, this function returns
    an unevaluated InverseMellinTransform object.

    Note that this function will assume x to be positive and real, regardless
    of the sympy assumptions!

    For a description of possible hints, refer to the docstring of
    :func:`sympy.transforms.IntegralTransform.doit`.

    >>> from sympy.integrals.transforms import inverse_mellin_transform
    >>> from sympy import oo, gamma
    >>> from sympy.abc import x, s
    >>> inverse_mellin_transform(gamma(s), s, x, (0, oo))
    exp(-x)

    The fundamental strip matters:

    >>> f = 1/(s**2 - 1)
    >>> inverse_mellin_transform(f, s, x, (-oo, -1))
    x*(1 - 1/x**2)*Heaviside(x - 1)/2
    >>> inverse_mellin_transform(f, s, x, (-1, 1))
    -x*Heaviside(-x + 1)/2 - Heaviside(x - 1)/(2*x)
    >>> inverse_mellin_transform(f, s, x, (1, oo))
    (-x**2/2 + 1/2)*Heaviside(-x + 1)/x

    See Also
    ========

    mellin_transform
    """
    return InverseMellinTransform(F, s, x, strip[0], strip[1]).doit(**hints)


##########################################################################
# Laplace Transform
##########################################################################

@_noconds
def _laplace_transform(f, t, s, simplify=True):
    """ The backend function for laplace transforms. """
    from sympy import (re, Max, exp, pi, Abs, Min, periodic_argument as arg,
                       cos, Wild, symbols)
    F = integrate(exp(-s*t) * f, (t, 0, oo))

    if not F.has(Integral):
        return _simplify(F, simplify), -oo, True

    if not F.is_Piecewise:
        raise IntegralTransformError('Laplace', f, 'could not compute integral')

    F, cond = F.args[0]
    if F.has(Integral):
        raise IntegralTransformError('Laplace', f, 'integral in unexpected form')

    a = -oo
    aux = True
    conds = conjuncts(to_cnf(cond))
    u = Dummy('u', real=True)
    p, q, w1, w2, w3 = symbols('p q w1 w2 w3', cls=Wild, exclude=[s])
    for c in conds:
        a_ = oo
        aux_ = []
        for d in disjuncts(c):
            m = d.match(abs(arg((s + w3)**p*q, w1)) < w2)
            if m:
                if m[q] > 0 and m[w2]/m[p] == pi/2:
                    d = re(s + m[w3]) > 0
            m = d.match(0 < cos(abs(arg(s, q)))*abs(s) - p)
            if m:
                d = re(s) > m[p]
            d_ = d.replace(re, lambda x: x.expand().as_real_imag()[0]).subs(re(s), t)
            if not d.is_Relational or (d.rel_op != '<' and d.rel_op != '<=') \
               or d_.has(s) or not d_.has(t):
                aux_ += [d]
                continue
            soln = _solve_inequality(d_, t)
            if not soln.is_Relational or \
               (soln.rel_op != '<' and soln.rel_op != '<='):
                aux_ += [d]
                continue
            if soln.lhs == t:
                raise IntegralTransformError('Laplace', f,
                                     'convergence not in half-plane?')
            else:
                a_ = Min(soln.lhs, a_)
        if a_ != oo:
            a = Max(a_, a)
        else:
            aux = And(aux, Or(*aux_))

    return _simplify(F, simplify), a, aux

class LaplaceTransform(IntegralTransform):
    """
    Class representing unevaluated laplace transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute laplace transforms, see the :func:`laplace_transform`
    docstring.
    """

    _name = 'Laplace'

    def _compute_transform(self, f, t, s, **hints):
        return _laplace_transform(f, t, s, **hints)

    def _as_integral(self, f, t, s):
        from sympy import Integral, exp
        return Integral(f*exp(-s*t), (t, 0, oo))

    """
    Class representing unevaluated laplace transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.
    For how to compute laplace transforms, see the :func:`laplace_transform`
    docstring.
    """
    def _collapse_extra(self, extra):
        from sympy import And, Max
        conds = []
        planes = []
        for plane, cond in extra:
            conds.append(cond)
            planes.append(plane)
        cond = And(*conds)
        plane = Max(*planes)
        if cond is False:
            raise IntegralTransformError('Laplace', None, 'No combined convergence.')
        return plane, cond

def laplace_transform(f, t, s, **hints):
    r"""
    Compute the Laplace Transform `F(s)` of `f(t)`,

    .. math :: F(s) = \int_0^\infty e^{-st} f(t) \mathrm{d}t.

    For all "sensible" functions, this converges absolutely in a
    half plane  `a < Re(s)`.

    This function returns (F, a, cond)
    where `F` is the laplace transform of `f`, `Re(s) > a` is the half-plane
    of convergence, and cond are auxiliary convergence conditions.

    If the integral cannot be computed in closed form, this function returns
    an unevaluated LaplaceTransform object.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.transforms.IntegralTransform.doit`. If ``noconds=True``,
    only `F` will be returned (i.e. not ``cond``, and also not the plane ``a``).

    >>> from sympy.integrals import laplace_transform
    >>> from sympy.abc import t, s, a
    >>> laplace_transform(t**a, t, s)[0:2]
    (s**(-a - 1)*gamma(a + 1), 0)

    See Also
    ========

    inverse_laplace_transform, mellin_transform, fourier_transform
    """
    return LaplaceTransform(f, t, s).doit(**hints)

@_noconds_(True)
def _inverse_laplace_transform(F, s, t_, plane, simplify=True):
    """ The backend function for inverse laplace transforms. """
    from sympy import exp, Heaviside, log, expand_complex, Integral, Piecewise
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
        return Heaviside(1/abs(coeff) - t**exponent)*e1 \
             + Heaviside(t**exponent - 1/abs(coeff))*e2

    try:
        f, cond = inverse_mellin_transform(F, s, exp(-t), (None, oo),
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
            cond = True
        f = f.replace(Piecewise, pw_simp)

    if f.is_Piecewise:
        # many of the functions called below can't work with piecewise
        # (b/c it has a bool in args)
        return f.subs(t, t_), cond

    u = Dummy('u')
    def simp_heaviside(arg):
        a = arg.subs(exp(-t), u)
        if a.has(t):
            return Heaviside(arg)
        rel = _solve_inequality(a > 0, u)
        if rel.lhs == u:
            k = log(rel.rhs)
            return Heaviside(t + k)
        else:
            k = log(rel.lhs)
            return Heaviside(-(t + k))
    f = f.replace(Heaviside, simp_heaviside)

    def simp_exp(arg):
        return expand_complex(exp(arg))
    f = f.replace(exp, simp_exp)

    # TODO it would be nice to fix cosh and sinh ... simplify messes these
    #      exponentials up

    return _simplify(f.subs(t, t_), simplify), cond

class InverseLaplaceTransform(IntegralTransform):
    """
    Class representing unevaluated inverse laplace transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute inverse laplace transforms, see the
    :func:`inverse_laplace_transform` docstring.
    """

    nargs = 4

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
        from sympy import I, Integral, exp
        c = self.__class__._c
        return Integral(exp(s*t)*F, (s, c - I*oo, c + I*oo))

def inverse_laplace_transform(F, s, t, plane=None, **hints):
    r"""
    Compute the inverse laplace transform of `F(s)`, defined as

    .. math :: f(t) = \int_{c-i\infty}^{c+i\infty} e^{st} F(s) \mathrm{d}s,

    for `c` so large that `F(s)` has no singularites in the
    half-plane `Re(s) > c-\epsilon`.

    The plane can be specified by
    argument ``plane``, but will be inferred if passed as None.

    Under certain regularity conditions, this recovers `f(t)` from its
    Laplace Transform `F(s)`, for non-negative `t`, and vice
    versa.

    If the integral cannot be computed in closed form, this function returns
    an unevaluated InverseLaplaceTransform object.

    Note that this function will always assume `t` to be real,
    regardless of the sympy assumption on `t`.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.transforms.IntegralTransform.doit`.

    >>> from sympy.integrals.transforms import inverse_laplace_transform
    >>> from sympy import exp, Symbol
    >>> from sympy.abc import s, t
    >>> a = Symbol('a', positive=True)
    >>> inverse_laplace_transform(exp(-a*s)/s, s, t)
    Heaviside(-a + t)

    See Also
    ========

    laplace_transform
    """
    return InverseLaplaceTransform(F, s, t, plane).doit(**hints)


##########################################################################
# Fourier Transform
##########################################################################

@_noconds_(True)
def _fourier_transform(f, x, k, a, b, name, simplify=True):
    """
    Compute a general fourier-type transform
        F(k) = a int_-oo^oo exp(b*I*x*k) f(x) dx.

    For suitable choice of a and b, this reduces to the standard fourier
    and inverse fourier transforms.
    """
    from sympy import exp, I, oo
    F = integrate(a*f*exp(b*I*x*k), (x, -oo, oo))

    if not F.has(Integral):
        return _simplify(F, simplify), True

    if not F.is_Piecewise:
        raise IntegralTransformError(name, f, 'could not compute integral')

    F, cond = F.args[0]
    if F.has(Integral):
        raise IntegralTransformError(name, f, 'integral in unexpected form')

    return _simplify(F, simplify), cond

class FourierTypeTransform(IntegralTransform):
    """ Base class for fourier transforms.
        Specify cls._a and cls._b.
    """

    def _compute_transform(self, f, x, k, **hints):
        return _fourier_transform(f, x, k,
                                  self.__class__._a, self.__class__._b,
                                  self.__class__._name, **hints)

    def _as_integral(self, f, x, k):
        from sympy import Integral, exp, I
        a = self.__class__._a
        b = self.__class__._b
        return Integral(a*f*exp(b*I*x*k), (x, -oo, oo))

class FourierTransform(FourierTypeTransform):
    """
    Class representing unevaluated fourier transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute fourier transforms, see the :func:`fourier_transform`
    docstring.
    """

    _name = 'Fourier'
    _a = 1
    _b = -2*S.Pi

def fourier_transform(f, x, k, **hints):
    r"""
    Compute the unitary, ordinary-frequency fourier transform of `f`, defined
    as

    .. math:: F(k) = \int_{-\infty}^\infty f(x) e^{-2\pi i x k} \mathrm{d} x.

    If the transform cannot be computed in closed form, this
    function returns an unevaluated FourierTransform object.

    For other fourier transform conventions, see the function
    :func:`sympy.integrals.transforms._fourier_transform`.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.transforms.IntegralTransform.doit`.
    Note that for this transform, by default ``noconds=True``.

    >>> from sympy import fourier_transform, exp
    >>> from sympy.abc import x, k
    >>> fourier_transform(exp(-x**2), x, k)
    sqrt(pi)*exp(-pi**2*k**2)
    >>> fourier_transform(exp(-x**2), x, k, noconds=False)
    (sqrt(pi)*exp(-pi**2*k**2), True)

    See Also
    ========

    inverse_fourier_transform, mellin_transform, laplace_transform
    """
    return FourierTransform(f, x, k).doit(**hints)

class InverseFourierTransform(FourierTypeTransform):
    """
    Class representing unevaluated inverse fourier transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute inverse fourier transforms, see the
    :func:`inverse_fourier_transform` docstring.
    """

    _name = 'Inverse Fourier'
    _a = 1
    _b = 2*S.Pi

def inverse_fourier_transform(F, k, x, **hints):
    r"""
    Compute the unitary, ordinary-frequency inverse fourier transform of `F`,
    defined as

    .. math:: f(x) = \int_{-\infty}^\infty F(k) e^{2\pi i x k} \mathrm{d} k.

    If the transform cannot be computed in closed form, this
    function returns an unevaluated InverseFourierTransform object.

    For other fourier transform conventions, see the function
    :func:`sympy.integrals.transforms._fourier_transform`.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.transforms.IntegralTransform.doit`.
    Note that for this transform, by default ``noconds=True``.

    >>> from sympy import inverse_fourier_transform, exp, sqrt, pi
    >>> from sympy.abc import x, k
    >>> inverse_fourier_transform(sqrt(pi)*exp(-(pi*k)**2), k, x)
    exp(-x**2)
    >>> inverse_fourier_transform(sqrt(pi)*exp(-(pi*k)**2), k, x, noconds=False)
    (exp(-x**2), True)

    See Also
    ========

    fourier_transform
    """
    return InverseFourierTransform(F, k, x).doit(**hints)
