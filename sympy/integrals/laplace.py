"""Laplace Transforms"""
from sympy.core import S, pi, I
from sympy.core.add import Add
from sympy.core.cache import cacheit
from sympy.core.function import (
    AppliedUndef, Derivative, expand, expand_complex, expand_mul, expand_trig,
    Lambda, WildFunction, diff)
from sympy.core.mul import Mul, prod
from sympy.core.relational import _canonical, Ge, Gt, Lt, Unequality, Eq
from sympy.core.sorting import ordered
from sympy.core.symbol import Dummy, symbols, Wild
from sympy.functions.elementary.complexes import (
    re, im, arg, Abs, polar_lift, periodic_argument)
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.hyperbolic import cosh, coth, sinh, asinh
from sympy.functions.elementary.miscellaneous import Max, Min, sqrt
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.trigonometric import cos, sin, atan
from sympy.functions.special.bessel import besseli, besselj, besselk, bessely
from sympy.functions.special.delta_functions import DiracDelta, Heaviside
from sympy.functions.special.error_functions import erf, erfc, Ei
from sympy.functions.special.gamma_functions import digamma, gamma, lowergamma
from sympy.integrals import integrate, Integral
from sympy.integrals.transforms import (
    _simplify, IntegralTransform, IntegralTransformError)
from sympy.logic.boolalg import to_cnf, conjuncts, disjuncts, Or, And
from sympy.matrices.matrices import MatrixBase
from sympy.polys.matrices.linsolve import _lin_eq2dict
from sympy.polys.polyerrors import PolynomialError
from sympy.polys.polyroots import roots
from sympy.polys.polytools import Poly
from sympy.polys.rationaltools import together
from sympy.polys.rootoftools import RootSum
from sympy.utilities.exceptions import (
    sympy_deprecation_warning, SymPyDeprecationWarning, ignore_warnings)
from sympy.utilities.misc import debug, debugf


def _simplifyconds(expr, s, a):
    r"""
    Naively simplify some conditions occurring in ``expr``,
    given that `\operatorname{Re}(s) > a`.

    Examples
    ========

    >>> from sympy.integrals.laplace import _simplifyconds
    >>> from sympy.abc import x
    >>> from sympy import sympify as S
    >>> _simplifyconds(abs(x**2) < 1, x, 1)
    False
    >>> _simplifyconds(abs(x**2) < 1, x, 2)
    False
    >>> _simplifyconds(abs(x**2) < 1, x, 0)
    Abs(x**2) < 1
    >>> _simplifyconds(abs(1/x**2) < 1, x, 1)
    True
    >>> _simplifyconds(S(1) < abs(x), x, 1)
    True
    >>> _simplifyconds(S(1) < abs(1/x), x, 1)
    False

    >>> from sympy import Ne
    >>> _simplifyconds(Ne(1, x**3), x, 1)
    True
    >>> _simplifyconds(Ne(1, x**3), x, 2)
    True
    >>> _simplifyconds(Ne(1, x**3), x, 0)
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
            if n > 0 and (Abs(ex1) <= Abs(a)**n) == S.true:
                return False
            if n < 0 and (Abs(ex1) >= Abs(a)**n) == S.true:
                return True
        except TypeError:
            pass

    def replie(x, y):
        """ simplify x < y """
        if (not (x.is_positive or isinstance(x, Abs))
                or not (y.is_positive or isinstance(y, Abs))):
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


def _laplace_transform_integration(f, t, s_, simplify=True):
    """ The backend function for doing Laplace transforms by integration.

    This backend assumes that the frontend has already split sums
    such that `f` is to an addition anymore.
    """
    s = Dummy('s')
    debugf('[LT _l_t_i ] started with (%s, %s, %s)', (f, t, s))
    debugf('[LT _l_t_i ]     and simplify=%s', (simplify, ))

    if f.has(DiracDelta):
        return None

    F = integrate(f*exp(-s*t), (t, S.Zero, S.Infinity))
    debugf('[LT _l_t_i ]     integrated: %s', (F, ))

    if not F.has(Integral):
        return _simplify(F.subs(s, s_), simplify), S.NegativeInfinity, S.true

    if not F.is_Piecewise:
        debug('[LT _l_t_i ]     not piecewise.')
        return None

    F, cond = F.args[0]
    if F.has(Integral):
        debug('[LT _l_t_i ]     integral in unexpected form.')
        return None

    def process_conds(conds):
        """ Turn ``conds`` into a strip and auxiliary conditions. """
        from sympy.solvers.inequalities import _solve_inequality
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
                if m and m[q].is_positive and m[w2]/m[p] == pi/2:
                    d = -re(s + m[w3]) < 0
                m = d.match(p - cos(w1*Abs(arg(s*w5))*w2)*Abs(s**w3)**w4 < 0)
                if not m:
                    m = d.match(
                        cos(p - Abs(periodic_argument(s**w1*w5, q))*w2) *
                        Abs(s**w3)**w4 < 0)
                if not m:
                    m = d.match(
                        p - cos(
                            Abs(periodic_argument(polar_lift(s)**w1*w5, q))*w2
                            )*Abs(s**w3)**w4 < 0)
                if m and all(m[wild].is_positive for wild in [
                        w1, w2, w3, w4, w5]):
                    d = re(s) > m[p]
                d_ = d.replace(
                    re, lambda x: x.expand().as_real_imag()[0]).subs(re(s), t)
                if (
                        not d.is_Relational or d.rel_op in ('==', '!=')
                        or d_.has(s) or not d_.has(t)):
                    aux_ += [d]
                    continue
                soln = _solve_inequality(d_, t)
                if not soln.is_Relational or soln.rel_op in ('==', '!='):
                    aux_ += [d]
                    continue
                if soln.lts == t:
                    debug('[LT _l_t_i ]     convergence not in half-plane.')
                    return None
                else:
                    a_ = Min(soln.lts, a_)
            if a_ is not S.Infinity:
                a = Max(a_, a)
            else:
                aux = And(aux, Or(*aux_))
        return a, aux.canonical if aux.is_Relational else aux

    conds = [process_conds(c) for c in disjuncts(cond)]
    conds2 = [x for x in conds if x[1] !=
              S.false and x[0] is not S.NegativeInfinity]
    if not conds2:
        conds2 = [x for x in conds if x[1] != S.false]
    conds = list(ordered(conds2))

    def cnt(expr):
        if expr in (True, False):
            return 0
        return expr.count_ops()
    conds.sort(key=lambda x: (-x[0], cnt(x[1])))

    if not conds:
        debug('[LT _l_t_i ]     no convergence found.')
        return None
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
        args = [_laplace_deep_collect(arg, t) for arg in args]
        if func.is_Add:
            return func(*args).collect(t)
        else:
            return func(*args)


@cacheit
def _laplace_build_rules():
    """
    This is an internal helper function that returns the table of Laplace
    transform rules in terms of the time variable `t` and the frequency
    variable `s`.  It is used by ``_laplace_apply_rules``.  Each entry is a
    tuple containing:

        (time domain pattern,
         frequency-domain replacement,
         condition for the rule to be applied,
         convergence plane,
         preparation function)

    The preparation function is a function with one argument that is applied
    to the expression before matching. For most rules it should be
    ``_laplace_deep_collect``.
    """
    t = Dummy('t')
    s = Dummy('s')
    a = Wild('a', exclude=[t])
    b = Wild('b', exclude=[t])
    n = Wild('n', exclude=[t])
    tau = Wild('tau', exclude=[t])
    omega = Wild('omega', exclude=[t])
    def dco(f): return _laplace_deep_collect(f, t)
    debug('_laplace_build_rules is building rules')

    laplace_transform_rules = [
        (a, a/s,
         S.true, S.Zero, dco),  # 4.2.1
        (DiracDelta(a*t-b), exp(-s*b/a)/Abs(a),
         Or(And(a > 0, b >= 0), And(a < 0, b <= 0)),
         S.NegativeInfinity, dco),  # Not in Bateman54
        (DiracDelta(a*t-b), S(0),
         Or(And(a < 0, b >= 0), And(a > 0, b <= 0)),
         S.NegativeInfinity, dco),  # Not in Bateman54
        (Heaviside(a*t-b), exp(-s*b/a)/s,
         And(a > 0, b > 0), S.Zero, dco),  # 4.4.1
        (Heaviside(a*t-b), (1-exp(-s*b/a))/s,
         And(a < 0, b < 0), S.Zero, dco),  # 4.4.1
        (Heaviside(a*t-b), 1/s,
         And(a > 0, b <= 0), S.Zero, dco),  # 4.4.1
        (Heaviside(a*t-b), 0,
         And(a < 0, b > 0), S.Zero, dco),  # 4.4.1
        (t, 1/s**2,
         S.true, S.Zero, dco),  # 4.2.3
        (1/(a*t+b), -exp(-b/a*s)*Ei(-b/a*s)/a,
         Abs(arg(b/a)) < pi, S.Zero, dco),  # 4.2.6
        (1/sqrt(a*t+b), sqrt(a*pi/s)*exp(b/a*s)*erfc(sqrt(b/a*s))/a,
         Abs(arg(b/a)) < pi, S.Zero, dco),  # 4.2.18
        ((a*t+b)**(-S(3)/2),
         2*b**(-S(1)/2)-2*(pi*s/a)**(S(1)/2)*exp(b/a*s) * erfc(sqrt(b/a*s))/a,
         Abs(arg(b/a)) < pi, S.Zero, dco),  # 4.2.20
        (sqrt(t)/(t+b), sqrt(pi/s)-pi*sqrt(b)*exp(b*s)*erfc(sqrt(b*s)),
         Abs(arg(b)) < pi, S.Zero, dco),  # 4.2.22
        (1/(a*sqrt(t) + t**(3/2)), pi*a**(S(1)/2)*exp(a*s)*erfc(sqrt(a*s)),
         S.true, S.Zero, dco),  # Not in Bateman54
        (t**n, gamma(n+1)/s**(n+1),
         n > -1, S.Zero, dco),  # 4.3.1
        ((a*t+b)**n, lowergamma(n+1, b/a*s)*exp(-b/a*s)/s**(n+1)/a,
         And(n > -1, Abs(arg(b/a)) < pi), S.Zero, dco),  # 4.3.4
        (t**n/(t+a), a**n*gamma(n+1)*lowergamma(-n, a*s),
         And(n > -1, Abs(arg(a)) < pi), S.Zero, dco),  # 4.3.7
        (exp(a*t-tau), exp(-tau)/(s-a),
         S.true, re(a), dco),  # 4.5.1
        (t*exp(a*t-tau), exp(-tau)/(s-a)**2,
         S.true, re(a), dco),  # 4.5.2
        (t**n*exp(a*t), gamma(n+1)/(s-a)**(n+1),
         re(n) > -1, re(a), dco),  # 4.5.3
        (exp(-a*t**2), sqrt(pi/4/a)*exp(s**2/4/a)*erfc(s/sqrt(4*a)),
         re(a) > 0, S.Zero, dco),  # 4.5.21
        (t*exp(-a*t**2),
         1/(2*a)-2/sqrt(pi)/(4*a)**(S(3)/2)*s*erfc(s/sqrt(4*a)),
         re(a) > 0, S.Zero, dco),  # 4.5.22
        (exp(-a/t), 2*sqrt(a/s)*besselk(1, 2*sqrt(a*s)),
         re(a) >= 0, S.Zero, dco),  # 4.5.25
        (sqrt(t)*exp(-a/t),
         S(1)/2*sqrt(pi/s**3)*(1+2*sqrt(a*s))*exp(-2*sqrt(a*s)),
         re(a) >= 0, S.Zero, dco),  # 4.5.26
        (exp(-a/t)/sqrt(t), sqrt(pi/s)*exp(-2*sqrt(a*s)),
         re(a) >= 0, S.Zero, dco),  # 4.5.27
        (exp(-a/t)/(t*sqrt(t)), sqrt(pi/a)*exp(-2*sqrt(a*s)),
         re(a) > 0, S.Zero, dco),  # 4.5.28
        (t**n*exp(-a/t), 2*(a/s)**((n+1)/2)*besselk(n+1, 2*sqrt(a*s)),
         re(a) > 0, S.Zero, dco),  # 4.5.29
        (exp(-2*sqrt(a*t)),
         s**(-1)-sqrt(pi*a)*s**(-S(3)/2)*exp(a/s) * erfc(sqrt(a/s)),
         Abs(arg(a)) < pi, S.Zero, dco),  # 4.5.31
        (exp(-2*sqrt(a*t))/sqrt(t), (pi/s)**(S(1)/2)*exp(a/s)*erfc(sqrt(a/s)),
         Abs(arg(a)) < pi, S.Zero, dco),  # 4.5.33
        (log(a*t), -log(exp(S.EulerGamma)*s/a)/s,
         a > 0, S.Zero, dco),  # 4.6.1
        (log(1+a*t), -exp(s/a)/s*Ei(-s/a),
         Abs(arg(a)) < pi, S.Zero, dco),  # 4.6.4
        (log(a*t+b), (log(b)-exp(s/b/a)/s*a*Ei(-s/b))/s*a,
         And(a > 0, Abs(arg(b)) < pi), S.Zero, dco),  # 4.6.5
        (log(t)/sqrt(t), -sqrt(pi/s)*log(4*s*exp(S.EulerGamma)),
         S.true, S.Zero, dco),  # 4.6.9
        (t**n*log(t), gamma(n+1)*s**(-n-1)*(digamma(n+1)-log(s)),
         re(n) > -1, S.Zero, dco),  # 4.6.11
        (log(a*t)**2, (log(exp(S.EulerGamma)*s/a)**2+pi**2/6)/s,
         a > 0, S.Zero, dco),  # 4.6.13
        (sin(omega*t), omega/(s**2+omega**2),
         S.true, Abs(im(omega)), dco),  # 4,7,1
        (Abs(sin(omega*t)), omega/(s**2+omega**2)*coth(pi*s/2/omega),
         omega > 0, S.Zero, dco),  # 4.7.2
        (sin(omega*t)/t, atan(omega/s),
         S.true, Abs(im(omega)), dco),  # 4.7.16
        (sin(omega*t)**2/t, log(1+4*omega**2/s**2)/4,
         S.true, 2*Abs(im(omega)), dco),  # 4.7.17
        (sin(omega*t)**2/t**2,
         omega*atan(2*omega/s)-s*log(1+4*omega**2/s**2)/4,
         S.true, 2*Abs(im(omega)), dco),  # 4.7.20
        (sin(2*sqrt(a*t)), sqrt(pi*a)/s/sqrt(s)*exp(-a/s),
         S.true, S.Zero, dco),  # 4.7.32
        (sin(2*sqrt(a*t))/t, pi*erf(sqrt(a/s)),
         S.true, S.Zero, dco),  # 4.7.34
        (cos(omega*t), s/(s**2+omega**2),
         S.true, Abs(im(omega)), dco),  # 4.7.43
        (cos(omega*t)**2, (s**2+2*omega**2)/(s**2+4*omega**2)/s,
         S.true, 2*Abs(im(omega)), dco),  # 4.7.45
        (sqrt(t)*cos(2*sqrt(a*t)), sqrt(pi)/2*s**(-S(5)/2)*(s-2*a)*exp(-a/s),
         S.true, S.Zero, dco),  # 4.7.66
        (cos(2*sqrt(a*t))/sqrt(t), sqrt(pi/s)*exp(-a/s),
         S.true, S.Zero, dco),  # 4.7.67
        (sin(a*t)*sin(b*t), 2*a*b*s/(s**2+(a+b)**2)/(s**2+(a-b)**2),
         S.true, Abs(im(a))+Abs(im(b)), dco),  # 4.7.78
        (cos(a*t)*sin(b*t), b*(s**2-a**2+b**2)/(s**2+(a+b)**2)/(s**2+(a-b)**2),
         S.true, Abs(im(a))+Abs(im(b)), dco),  # 4.7.79
        (cos(a*t)*cos(b*t), s*(s**2+a**2+b**2)/(s**2+(a+b)**2)/(s**2+(a-b)**2),
         S.true, Abs(im(a))+Abs(im(b)), dco),  # 4.7.80
        (sinh(a*t), a/(s**2-a**2),
         S.true, Abs(re(a)), dco),  # 4.9.1
        (cosh(a*t), s/(s**2-a**2),
         S.true, Abs(re(a)), dco),  # 4.9.2
        (sinh(a*t)**2, 2*a**2/(s**3-4*a**2*s),
         S.true, 2*Abs(re(a)), dco),  # 4.9.3
        (cosh(a*t)**2, (s**2-2*a**2)/(s**3-4*a**2*s),
         S.true, 2*Abs(re(a)), dco),  # 4.9.4
        (sinh(a*t)/t, log((s+a)/(s-a))/2,
         S.true, Abs(re(a)), dco),  # 4.9.12
        (t**n*sinh(a*t), gamma(n+1)/2*((s-a)**(-n-1)-(s+a)**(-n-1)),
         n > -2, Abs(a), dco),  # 4.9.18
        (t**n*cosh(a*t), gamma(n+1)/2*((s-a)**(-n-1)+(s+a)**(-n-1)),
         n > -1, Abs(a), dco),  # 4.9.19
        (sinh(2*sqrt(a*t)), sqrt(pi*a)/s/sqrt(s)*exp(a/s),
         S.true, S.Zero, dco),  # 4.9.34
        (cosh(2*sqrt(a*t)), 1/s+sqrt(pi*a)/s/sqrt(s)*exp(a/s)*erf(sqrt(a/s)),
         S.true, S.Zero, dco),  # 4.9.35
        (
            sqrt(t)*sinh(2*sqrt(a*t)),
            pi**(S(1)/2)*s**(-S(5)/2)*(s/2+a) *
            exp(a/s)*erf(sqrt(a/s))-a**(S(1)/2)*s**(-2),
            S.true, S.Zero, dco),  # 4.9.36
        (sqrt(t)*cosh(2*sqrt(a*t)), pi**(S(1)/2)*s**(-S(5)/2)*(s/2+a)*exp(a/s),
         S.true, S.Zero, dco),  # 4.9.37
        (sinh(2*sqrt(a*t))/sqrt(t),
         pi**(S(1)/2)*s**(-S(1)/2)*exp(a/s) * erf(sqrt(a/s)),
            S.true, S.Zero, dco),  # 4.9.38
        (cosh(2*sqrt(a*t))/sqrt(t), pi**(S(1)/2)*s**(-S(1)/2)*exp(a/s),
         S.true, S.Zero, dco),  # 4.9.39
        (sinh(sqrt(a*t))**2/sqrt(t), pi**(S(1)/2)/2*s**(-S(1)/2)*(exp(a/s)-1),
         S.true, S.Zero, dco),  # 4.9.40
        (cosh(sqrt(a*t))**2/sqrt(t), pi**(S(1)/2)/2*s**(-S(1)/2)*(exp(a/s)+1),
         S.true, S.Zero, dco),  # 4.9.41
        (erf(a*t), exp(s**2/(2*a)**2)*erfc(s/(2*a))/s,
         4*Abs(arg(a)) < pi, S.Zero, dco),  # 4.12.2
        (erf(sqrt(a*t)), sqrt(a)/sqrt(s+a)/s,
         S.true, Max(S.Zero, -re(a)), dco),  # 4.12.4
        (exp(a*t)*erf(sqrt(a*t)), sqrt(a)/sqrt(s)/(s-a),
         S.true, Max(S.Zero, re(a)), dco),  # 4.12.5
        (erf(sqrt(a/t)/2), (1-exp(-sqrt(a*s)))/s,
         re(a) > 0, S.Zero, dco),  # 4.12.6
        (erfc(sqrt(a*t)), (sqrt(s+a)-sqrt(a))/sqrt(s+a)/s,
         S.true, -re(a), dco),  # 4.12.9
        (exp(a*t)*erfc(sqrt(a*t)), 1/(s+sqrt(a*s)),
         S.true, S.Zero, dco),  # 4.12.10
        (erfc(sqrt(a/t)/2), exp(-sqrt(a*s))/s,
         re(a) > 0, S.Zero, dco),  # 4.2.11
        (besselj(n, a*t), a**n/(sqrt(s**2+a**2)*(s+sqrt(s**2+a**2))**n),
         re(n) > -1, Abs(im(a)), dco),  # 4.14.1
        (t**b*besselj(n, a*t),
         2**n/sqrt(pi)*gamma(n+S.Half)*a**n*(s**2+a**2)**(-n-S.Half),
         And(re(n) > -S.Half, Eq(b, n)), Abs(im(a)), dco),  # 4.14.7
        (t**b*besselj(n, a*t),
         2**(n+1)/sqrt(pi)*gamma(n+S(3)/2)*a**n*s*(s**2+a**2)**(-n-S(3)/2),
         And(re(n) > -1, Eq(b, n+1)), Abs(im(a)), dco),  # 4.14.8
        (besselj(0, 2*sqrt(a*t)), exp(-a/s)/s,
         S.true, S.Zero, dco),  # 4.14.25
        (t**(b)*besselj(n, 2*sqrt(a*t)), a**(n/2)*s**(-n-1)*exp(-a/s),
         And(re(n) > -1, Eq(b, n*S.Half)), S.Zero, dco),  # 4.14.30
        (besselj(0, a*sqrt(t**2+b*t)),
         exp(b*s-b*sqrt(s**2+a**2))/sqrt(s**2+a**2),
         Abs(arg(b)) < pi, Abs(im(a)), dco),  # 4.15.19
        (besseli(n, a*t), a**n/(sqrt(s**2-a**2)*(s+sqrt(s**2-a**2))**n),
         re(n) > -1, Abs(re(a)), dco),  # 4.16.1
        (t**b*besseli(n, a*t),
         2**n/sqrt(pi)*gamma(n+S.Half)*a**n*(s**2-a**2)**(-n-S.Half),
         And(re(n) > -S.Half, Eq(b, n)), Abs(re(a)), dco),  # 4.16.6
        (t**b*besseli(n, a*t),
         2**(n+1)/sqrt(pi)*gamma(n+S(3)/2)*a**n*s*(s**2-a**2)**(-n-S(3)/2),
         And(re(n) > -1, Eq(b, n+1)), Abs(re(a)), dco),  # 4.16.7
        (t**(b)*besseli(n, 2*sqrt(a*t)), a**(n/2)*s**(-n-1)*exp(a/s),
         And(re(n) > -1, Eq(b, n*S.Half)), S.Zero, dco),  # 4.16.18
        (bessely(0, a*t), -2/pi*asinh(s/a)/sqrt(s**2+a**2),
         S.true, Abs(im(a)), dco),  # 4.15.44
        (besselk(0, a*t), log((s + sqrt(s**2-a**2))/a)/(sqrt(s**2-a**2)),
         S.true, -re(a), dco)  # 4.16.23
    ]
    return laplace_transform_rules, t, s


def _laplace_rule_timescale(f, t, s):
    """
    This function applies the time-scaling rule of the Laplace transform in
    a straight-forward way. For example, if it gets ``(f(a*t), t, s)``, it will
    compute ``LaplaceTransform(f(t)/a, t, s/a)`` if ``a>0``.
    """

    a = Wild('a', exclude=[t])
    g = WildFunction('g', nargs=1)
    ma1 = f.match(g)
    if ma1:
        arg = ma1[g].args[0].collect(t)
        ma2 = arg.match(a*t)
        if ma2 and ma2[a].is_positive and ma2[a] != 1:
            debug('_laplace_apply_prog rules match:')
            debugf('      f:    %s _ %s, %s )', (f, ma1, ma2))
            debug('      rule: time scaling (4.1.4)')
            r, pr, cr = _laplace_transform(1/ma2[a]*ma1[g].func(t),
                                           t, s/ma2[a], simplify=False)
            return (r, pr, cr)
    return None


def _laplace_rule_heaviside(f, t, s):
    """
    This function deals with time-shifted Heaviside step functions. If the time
    shift is positive, it applies the time-shift rule of the Laplace transform.
    For example, if it gets ``(Heaviside(t-a)*f(t), t, s)``, it will compute
    ``exp(-a*s)*LaplaceTransform(f(t+a), t, s)``.

    If the time shift is negative, the Heaviside function is simply removed
    as it means nothing to the Laplace transform.

    The function does not remove a factor ``Heaviside(t)``; this is done by
    the simple rules.
    """

    a = Wild('a', exclude=[t])
    y = Wild('y')
    g = Wild('g')
    ma1 = f.match(Heaviside(y)*g)
    if ma1:
        ma2 = ma1[y].match(t-a)
        if ma2 and ma2[a].is_positive:
            debug('_laplace_apply_prog_rules match:')
            debugf('      f:    %s ( %s, %s )', (f, ma1, ma2))
            debug('      rule: time shift (4.1.4)')
            r, pr, cr = _laplace_transform(ma1[g].subs(t, t+ma2[a]), t, s,
                                           simplify=False)
            return (exp(-ma2[a]*s)*r, pr, cr)
        if ma2 and ma2[a].is_negative:
            debug('_laplace_apply_prog_rules match:')
            debugf('      f:    %s ( %s, %s )', (f, ma1, ma2))
            debug('      rule: Heaviside factor, negative time shift (4.1.4)')
            r, pr, cr = _laplace_transform(ma1[g], t, s, simplify=False)
            return (r, pr, cr)
    return None


def _laplace_rule_exp(f, t, s):
    """
    If this function finds a factor ``exp(a*t)``, it applies the
    frequency-shift rule of the Laplace transform and adjusts the convergence
    plane accordingly.  For example, if it gets ``(exp(-a*t)*f(t), t, s)``, it
    will compute ``LaplaceTransform(f(t), t, s+a)``.
    """

    a = Wild('a', exclude=[t])
    y = Wild('y')
    z = Wild('z')
    ma1 = f.match(exp(y)*z)
    if ma1:
        ma2 = ma1[y].collect(t).match(a*t)
        if ma2:
            debug('_laplace_apply_prog_rules match:')
            debugf('      f:    %s ( %s, %s )', (f, ma1, ma2))
            debug('      rule: multiply with exp (4.1.5)')
            r, pr, cr = _laplace_transform(ma1[z], t, s-ma2[a],
                                           simplify=False)
            return (r, pr+re(ma2[a]), cr)
    return None


def _laplace_rule_delta(f, t, s):
    """
    If this function finds a factor ``DiracDelta(b*t-a)``, it applies the
    masking property of the delta distribution. For example, if it gets
    ``(DiracDelta(t-a)*f(t), t, s)``, it will return
    ``(f(a)*exp(-a*s), -a, True)``.
    """
    # This rule is not in Bateman54

    a = Wild('a', exclude=[t])
    b = Wild('b', exclude=[t])

    y = Wild('y')
    z = Wild('z')
    ma1 = f.match(DiracDelta(y)*z)
    if ma1 and not ma1[z].has(DiracDelta):
        ma2 = ma1[y].collect(t).match(b*t-a)
        if ma2:
            debug('_laplace_apply_prog_rules match:')
            debugf('      f:    %s ( %s, %s )', (f, ma1, ma2))
            debug('      rule: multiply with DiracDelta')
            loc = ma2[a]/ma2[b]
            if re(loc) >= 0 and im(loc) == 0:
                r = exp(-ma2[a]/ma2[b]*s)*ma1[z].subs(t, ma2[a]/ma2[b])/ma2[b]
                return (r, S.NegativeInfinity, S.true)
            else:
                return (0, S.NegativeInfinity, S.true)
        if ma1[y].is_polynomial(t):
            ro = roots(ma1[y], t)
            if roots is not {} and set(ro.values()) == {1}:
                slope = diff(ma1[y], t)
                r = Add(
                    *[exp(-x*s)*ma1[z].subs(t, s)/slope.subs(t, x)
                      for x in list(ro.keys()) if im(x) == 0 and re(x) >= 0])
                return (r, S.NegativeInfinity, S.true)
    return None


def _laplace_trig_split(fn):
    """
    Helper function for `_laplace_rule_trig`.  This function returns two terms
    `f` and `g`.  `f` contains all product terms with sin, cos, sinh, cosh in
    them; `g` contains everything else.
    """
    trigs = [S.One]
    other = [S.One]
    for term in Mul.make_args(fn):
        if term.has(sin, cos, sinh, cosh, exp):
            trigs.append(term)
        else:
            other.append(term)
    f = Mul(*trigs)
    g = Mul(*other)
    return f, g


def _laplace_trig_expsum(f, t):
    """
    Helper function for `_laplace_rule_trig`.  This function expects the `f`
    from `_laplace_trig_split`.  It returns two lists `xm` and `xn`.  `xm` is
    a list of dictionaries with keys `k` and `a` representing a function
    `k*exp(a*t)`.  `xn` is a list of all terms that cannot be brought into
    that form, which may happen, e.g., when a trigonometric function has
    another function in its argument.
    """
    m = Wild('m')
    p = Wild('p', exclude=[t])
    xm = []
    xn = []

    x1 = f.rewrite(exp).expand()

    for term in Add.make_args(x1):
        if not term.has(t):
            xm.append({'k': term, 'a': 0, re: 0, im: 0})
            continue
        term = term.powsimp(combine='exp')
        if (r := term.match(p*exp(m))) is not None:
            if (mp := r[m].as_poly(t)) is not None:
                mc = mp.all_coeffs()
                if len(mc) == 2:
                    xm.append({
                        'k': r[p]*exp(mc[1]), 'a': mc[0],
                        re: re(mc[0]), im: im(mc[0])})
                else:
                    xn.append(term)
            else:
                xn.append(term)
        else:
            xn.append(term)
    return xm, xn


def _laplace_trig_ltex(xm, t, s):
    """
    Helper function for `_laplace_rule_trig`.  This function takes the list of
    exponentials `xm` from `_laplace_trig_expsum` and simplifies complex
    conjugate and real symmetric poles.  It returns the result as a sum and
    the convergence plane.
    """
    results = []
    planes = []

    def _simpc(coeffs):
        nc = coeffs.copy()
        for k in range(len(nc)):
            ri = nc[k].as_real_imag()
            if ri[0].has(im):
                nc[k] = nc[k].rewrite(cos)
            else:
                nc[k] = (ri[0] + I*ri[1]).rewrite(cos)
        return nc

    def _quadpole(t1, k1, k2, k3, s):
        a, k0, a_r, a_i = t1['a'], t1['k'], t1[re], t1[im]
        nc = [
            k0 + k1 + k2 + k3,
            a*(k0 + k1 - k2 - k3) - 2*I*a_i*k1 + 2*I*a_i*k2,
            (
                a**2*(-k0 - k1 - k2 - k3) +
                a*(4*I*a_i*k0 + 4*I*a_i*k3) +
                4*a_i**2*k0 + 4*a_i**2*k3),
            (
                a**3*(-k0 - k1 + k2 + k3) +
                a**2*(4*I*a_i*k0 + 2*I*a_i*k1 - 2*I*a_i*k2 - 4*I*a_i*k3) +
                a*(4*a_i**2*k0 - 4*a_i**2*k3))
        ]
        dc = [
            S.One, S.Zero, 2*a_i**2 - 2*a_r**2,
            S.Zero, a_i**4 + 2*a_i**2*a_r**2 + a_r**4]
        n = Add(
            *[x*s**y for x, y in zip(_simpc(nc), range(len(nc))[::-1])])
        d = Add(
            *[x*s**y for x, y in zip(dc, range(len(dc))[::-1])])
        debugf('        quadpole: (%s) / (%s)', (n, d))
        return n/d

    def _ccpole(t1, k1, s):
        a, k0, a_r, a_i = t1['a'], t1['k'], t1[re], t1[im]
        nc = [k0 + k1, -a*k0 - a*k1 + 2*I*a_i*k0]
        dc = [S.One, -2*a_r, a_i**2 + a_r**2]
        n = Add(
            *[x*s**y for x, y in zip(_simpc(nc), range(len(nc))[::-1])])
        d = Add(
            *[x*s**y for x, y in zip(dc, range(len(dc))[::-1])])
        debugf('        ccpole: (%s) / (%s)', (n, d))
        return n/d

    def _rspole(t1, k2, s):
        a, k0, a_r, a_i = t1['a'], t1['k'], t1[re], t1[im]
        nc = [k0 + k2, a*k0 - a*k2 - 2*I*a_i*k0]
        dc = [S.One, -2*I*a_i, -a_i**2 - a_r**2]
        n = Add(
            *[x*s**y for x, y in zip(_simpc(nc), range(len(nc))[::-1])])
        d = Add(
            *[x*s**y for x, y in zip(dc, range(len(dc))[::-1])])
        debugf('        rspole: (%s) / (%s)', (n, d))
        return n/d

    def _sypole(t1, k3, s):
        a, k0 = t1['a'], t1['k']
        nc = [k0 + k3, a*(k0 - k3)]
        dc = [S.One, S.Zero, -a**2]
        n = Add(
            *[x*s**y for x, y in zip(_simpc(nc), range(len(nc))[::-1])])
        d = Add(
            *[x*s**y for x, y in zip(dc, range(len(dc))[::-1])])
        debugf('        sypole: (%s) / (%s)', (n, d))
        return n/d

    def _simplepole(t1, s):
        a, k0 = t1['a'], t1['k']
        n = k0
        d = s - a
        debugf('        simplepole: (%s) / (%s)', (n, d))
        return n/d

    while len(xm) > 0:
        t1 = xm.pop()
        i_imagsym = None
        i_realsym = None
        i_pointsym = None
        # The following code checks all remaining poles. If t1 is a pole at
        # a+b*I, then we check for a-b*I, -a+b*I, and -a-b*I, and
        # assign the respective indices to i_imagsym, i_realsym, i_pointsym.
        # -a-b*I / i_pointsym only applies if both a and b are != 0.
        for i in range(len(xm)):
            real_eq = t1[re] == xm[i][re]
            realsym = t1[re] == -xm[i][re]
            imag_eq = t1[im] == xm[i][im]
            imagsym = t1[im] == -xm[i][im]
            if realsym and imagsym and t1[re] != 0 and t1[im] != 0:
                i_pointsym = i
            elif realsym and imag_eq and t1[re] != 0:
                i_realsym = i
            elif real_eq and imagsym and t1[im] != 0:
                i_imagsym = i

        # The next part looks for four possible pole constellations:
        # quad:   a+b*I, a-b*I, -a+b*I, -a-b*I
        # cc:     a+b*I, a-b*I (a may be zero)
        # quad:   a+b*I, -a+b*I (b may be zero)
        # point:  a+b*I, -a-b*I (a!=0 and b!=0 is needed, but that has been
        #                        asserted when finding i_pointsym above.)
        # If none apply, then t1 is a simple pole.
        if (
                i_imagsym is not None and i_realsym is not None
                and i_pointsym is not None):
            results.append(
                _quadpole(t1,
                          xm[i_imagsym]['k'], xm[i_realsym]['k'],
                          xm[i_pointsym]['k'], s))
            planes.append(Abs(re(t1['a'])))
            # The three additional poles have now been used; to pop them
            # easily we have to do it from the back.
            indices_to_pop = [i_imagsym, i_realsym, i_pointsym]
            indices_to_pop.sort(reverse=True)
            for i in indices_to_pop:
                xm.pop(i)
        elif i_imagsym is not None:
            results.append(_ccpole(t1, xm[i_imagsym]['k'], s))
            planes.append(t1[re])
            xm.pop(i_imagsym)
        elif i_realsym is not None:
            results.append(_rspole(t1, xm[i_realsym]['k'], s))
            planes.append(Abs(t1[re]))
            xm.pop(i_realsym)
        elif i_pointsym is not None:
            results.append(_sypole(t1, xm[i_pointsym]['k'], s))
            planes.append(Abs(t1[re]))
            xm.pop(i_pointsym)
        else:
            results.append(_simplepole(t1, s))
            planes.append(t1[re])

    return Add(*results), Max(*planes)


def _laplace_rule_trig(fn, t_, s, doit=True, **hints):
    """
    This rule covers trigonometric factors by splitting everything into a
    sum of exponential functions and collecting complex conjugate poles and
    real symmetric poles.
    """
    t = Dummy('t', real=True)

    if not fn.has(sin, cos, sinh, cosh):
        return None

    debugf('_laplace_rule_trig: (%s, %s, %s)', (fn, t_, s))

    f, g = _laplace_trig_split(fn.subs(t_, t))
    debugf('    f = %s\n    g = %s', (f, g))

    xm, xn = _laplace_trig_expsum(f, t)
    debugf('    xm = %s\n    xn = %s', (xm, xn))

    if len(xn) > 0:
        # not implemented yet
        debug('    --> xn is not empty; giving up.')
        return None

    if not g.has(t):
        r, p = _laplace_trig_ltex(xm, t, s)
        return g*r, p, S.true
    else:
        # Just transform `g` and make frequency-shifted copies
        planes = []
        results = []
        G, G_plane, G_cond = _laplace_transform(g, t, s)
        for x1 in xm:
            results.append(x1['k']*G.subs(s, s-x1['a']))
            planes.append(G_plane+re(x1['a']))
    return Add(*results).subs(t, t_), Max(*planes), G_cond


def _laplace_rule_diff(f, t, s, doit=True, **hints):
    """
    This function looks for derivatives in the time domain and replaces it
    by factors of `s` and initial conditions in the frequency domain. For
    example, if it gets ``(diff(f(t), t), t, s)``, it will compute
    ``s*LaplaceTransform(f(t), t, s) - f(0)``.
    """

    a = Wild('a', exclude=[t])
    n = Wild('n', exclude=[t])
    g = WildFunction('g')
    ma1 = f.match(a*Derivative(g, (t, n)))
    if ma1 and ma1[n].is_integer:
        m = [z.has(t) for z in ma1[g].args]
        if sum(m) == 1:
            debug('_laplace_apply_rules match:')
            debugf('      f, n: %s, %s', (f, ma1[n]))
            debug('      rule: time derivative (4.1.8)')
            d = []
            for k in range(ma1[n]):
                if k == 0:
                    y = ma1[g].subs(t, 0)
                else:
                    y = Derivative(ma1[g], (t, k)).subs(t, 0)
                d.append(s**(ma1[n]-k-1)*y)
            r, pr, cr = _laplace_transform(ma1[g], t, s, simplify=False)
            return (ma1[a]*(s**ma1[n]*r - Add(*d)),  pr, cr)
    return None


def _laplace_rule_sdiff(f, t, s, doit=True, **hints):
    """
    This function looks for multiplications with polynoimials in `t` as they
    correspond to differentiation in the frequency domain. For example, if it
    gets ``(t*f(t), t, s)``, it will compute
    ``-Derivative(LaplaceTransform(f(t), t, s), s)``.
    """

    if f.is_Mul:
        pfac = [1]
        ofac = [1]
        for fac in Mul.make_args(f):
            if fac.is_polynomial(t):
                pfac.append(fac)
            else:
                ofac.append(fac)
        if len(pfac) > 1:
            pex = prod(pfac)
            pc = Poly(pex, t).all_coeffs()
            N = len(pc)
            if N > 1:
                debug('_laplace_apply_rules match:')
                debugf('      f, n: %s, %s', (f, pfac))
                debug('      rule: frequency derivative (4.1.6)')
                oex = prod(ofac)
                r_, p_, c_ = _laplace_transform(oex, t, s, simplify=False)
                deri = [r_]
                d1 = False
                try:
                    d1 = -diff(deri[-1], s)
                except ValueError:
                    d1 = False
                if r_.has(LaplaceTransform):
                    for k in range(N-1):
                        deri.append((-1)**(k+1)*Derivative(r_, s, k+1))
                else:
                    if d1:
                        deri.append(d1)
                        for k in range(N-2):
                            deri.append(-diff(deri[-1], s))
                if d1:
                    r = Add(*[pc[N-n-1]*deri[n] for n in range(N)])
                    return (r, p_, c_)
    return None


def _laplace_expand(f, t, s, doit=True, **hints):
    """
    This function tries to expand its argument with successively stronger
    methods: first it will expand on the top level, then it will expand any
    multiplications in depth, then it will try all avilable expansion methods,
    and finally it will try to expand trigonometric functions.

    If it can expand, it will then compute the Laplace transform of the
    expanded term.
    """

    if f.is_Add:
        return None
    r = expand(f, deep=False)
    if r.is_Add:
        return _laplace_transform(r, t, s, simplify=False)
    r = expand_mul(f)
    if r.is_Add:
        return _laplace_transform(r, t, s, simplify=False)
    r = expand(f)
    if r.is_Add:
        return _laplace_transform(r, t, s, simplify=False)
    if r != f:
        return _laplace_transform(r, t, s, simplify=False)
    r = expand(expand_trig(f))
    if r.is_Add:
        return _laplace_transform(r, t, s, simplify=False)
    return None


def _laplace_apply_prog_rules(f, t, s):
    """
    This function applies all program rules and returns the result if one
    of them gives a result.
    """

    prog_rules = [_laplace_rule_heaviside, _laplace_rule_delta,
                  _laplace_rule_timescale, _laplace_rule_exp,
                  _laplace_rule_trig,
                  _laplace_rule_diff, _laplace_rule_sdiff]

    for p_rule in prog_rules:
        if (L := p_rule(f, t, s)) is not None:
            return L
    return None


def _laplace_apply_simple_rules(f, t, s):
    """
    This function applies all simple rules and returns the result if one
    of them gives a result.
    """
    simple_rules, t_, s_ = _laplace_build_rules()
    prep_old = ''
    prep_f = ''
    for t_dom, s_dom, check, plane, prep in simple_rules:
        if prep_old != prep:
            prep_f = prep(f.subs({t: t_}))
            prep_old = prep
        ma = prep_f.match(t_dom)
        if ma:
            try:
                c = check.xreplace(ma)
            except TypeError:
                # This may happen if the time function has imaginary
                # numbers in it. Then we give up.
                continue
            if c == S.true:
                debug('_laplace_apply_simple_rules match:')
                debugf('      f:     %s', (f,))
                debugf('      rule:  %s o---o %s', (t_dom, s_dom))
                debugf('      match: %s', (ma, ))
                return (s_dom.xreplace(ma).subs({s_: s}),
                        plane.xreplace(ma), S.true)
    return None


def _laplace_transform(fn, t_, s_, simplify=True):
    """
    Front-end function of the Laplace transform. It tries to apply all known
    rules recursively, and if everything else fails, it tries to integrate.
    """
    debugf('[LT _l_t] (%s, %s, %s)', (fn, t_, s_))

    terms = Add.make_args(fn)
    terms_s = []
    planes = []
    conditions = []
    for ff in terms:
        k, ft = ff.as_independent(t_, as_Add=False)
        if (r := _laplace_apply_simple_rules(ft, t_, s_)) is not None:
            pass
        elif (r := _laplace_apply_prog_rules(ft, t_, s_)) is not None:
            pass
        elif (r := _laplace_expand(ft, t_, s_)) is not None:
            pass
        elif any(undef.has(t_) for undef in ft.atoms(AppliedUndef)):
            # If there are undefined functions f(t) then integration is
            # unlikely to do anything useful so we skip it and given an
            # unevaluated LaplaceTransform.
            r = (LaplaceTransform(ft, t_, s_), S.NegativeInfinity, True)
        elif (r := _laplace_transform_integration(
                ft, t_, s_, simplify=simplify)) is not None:
            pass
        else:
            r = (LaplaceTransform(ft, t_, s_), S.NegativeInfinity, True)
        (ri_, pi_, ci_) = r
        terms_s.append(k*ri_)
        planes.append(pi_)
        conditions.append(ci_)

    result = Add(*terms_s)
    if simplify:
        result = result.simplify(doit=False)
    plane = Max(*planes)
    condition = And(*conditions)

    return result, plane, condition


class LaplaceTransform(IntegralTransform):
    """
    Class representing unevaluated Laplace transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute Laplace transforms, see the :func:`laplace_transform`
    docstring.

    If this is called with ``.doit()``, it returns the Laplace transform as an
    expression. If it is called with ``.doit(noconds=False)``, it returns a
    tuple containing the same expression, a convergence plane, and conditions.
    """

    _name = 'Laplace'

    def _compute_transform(self, f, t, s, **hints):
        _simplify = hints.get('simplify', False)
        LT = _laplace_transform_integration(f, t, s, simplify=_simplify)
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
        if cond == S.false:
            raise IntegralTransformError(
                'Laplace', None, 'No combined convergence.')
        return plane, cond

    def doit(self, **hints):
        """
        Try to evaluate the transform in closed form.

        Explanation
        ===========

        Standard hints are the following:
        - ``noconds``:  if True, do not return convergence conditions. The
        default setting is `True`.
        - ``simplify``: if True, it simplifies the final result. The
        default setting is `False`.
        """
        _noconds = hints.get('noconds', True)
        _simplify = hints.get('simplify', False)

        debugf('[LT doit] (%s, %s, %s)', (self.function,
                                          self.function_variable,
                                          self.transform_variable))

        t_ = self.function_variable
        s_ = self.transform_variable
        fn = self.function

        r = _laplace_transform(fn, t_, s_, simplify=_simplify)

        if _noconds:
            return r[0]
        else:
            return r


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
    rules are applied, and whether integration is attempted, you can switch
    debug information on by setting ``sympy.SYMPY_DEBUG=True``. The numbers
    of the rules in the debug information (and the code) refer to Bateman's
    Tables of Integral Transforms [1].

    The lower bound is `0-`, meaning that this bound should be approached
    from the lower side. This is only necessary if distributions are involved.
    At present, it is only done if `f(t)` contains ``DiracDelta``, in which
    case the Laplace transform is computed implicitly as

    .. math ::
        F(s) = \lim_{\tau\to 0^{-}} \int_{\tau}^\infty e^{-st}
        f(t) \mathrm{d}t

    by applying rules.

    If the Laplace transform cannot be fully computed in closed form, this
    function returns expressions containing unevaluated
    :class:`LaplaceTransform` objects.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`. If
    ``noconds=True``, only `F` will be returned (i.e. not ``cond``, and also
    not the plane ``a``).

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
    >>> laplace_transform(DiracDelta(t)-a*exp(-a*t), t, s, simplify=True)
    (s/(a + s), -re(a), True)

    References
    ==========

    .. [1] Erdelyi, A. (ed.), Tables of Integral Transforms, Volume 1,
           Bateman Manuscript Prooject, McGraw-Hill (1954), available:
           https://resolver.caltech.edu/CaltechAUTHORS:20140123-101456353

    See Also
    ========

    inverse_laplace_transform, mellin_transform, fourier_transform
    hankel_transform, inverse_hankel_transform

    """

    _noconds = hints.get('noconds', False)
    _simplify = hints.get('simplify', False)

    if isinstance(f, MatrixBase) and hasattr(f, 'applyfunc'):

        conds = not hints.get('noconds', False)

        if conds and legacy_matrix:
            adt = 'deprecated-laplace-transform-matrix'
            sympy_deprecation_warning(
                """
Calling laplace_transform() on a Matrix with noconds=False (the default) is
deprecated. Either noconds=True or use legacy_matrix=False to get the new
behavior.
                """,
                deprecated_since_version='1.9',
                active_deprecations_target=adt,
            )
            # Temporarily disable the deprecation warning for non-Expr objects
            # in Matrix
            with ignore_warnings(SymPyDeprecationWarning):
                return f.applyfunc(
                    lambda fij: laplace_transform(fij, t, s, **hints))
        else:
            elements_trans = [laplace_transform(
                fij, t, s, **hints) for fij in f]
            if conds:
                elements, avals, conditions = zip(*elements_trans)
                f_laplace = type(f)(*f.shape, elements)
                return f_laplace, Max(*avals), And(*conditions)
            else:
                return type(f)(*f.shape, elements_trans)

    LT = LaplaceTransform(f, t, s).doit(noconds=False, simplify=_simplify)

    if not _noconds:
        return LT
    else:
        return LT[0]


def _inverse_laplace_transform_integration(F, s, t_, plane, simplify=True):
    """ The backend function for inverse Laplace transforms. """
    from sympy.integrals.meijerint import meijerint_inversion, _get_coeff_exp
    from sympy.integrals.transforms import inverse_mellin_transform

    # There are two strategies we can try:
    # 1) Use inverse mellin transform, related by a simple change of variables.
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
        return (
            Heaviside(1/Abs(coeff) - t**exponent)*e1 +
            Heaviside(t**exponent - 1/Abs(coeff))*e2)

    if F.is_rational_function(s):
        F = F.apart(s)

    if F.is_Add:
        f = Add(
            *[_inverse_laplace_transform_integration(X, s, t, plane, simplify)
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
            return None
        if f.is_Piecewise:
            f, cond = f.args[0]
            if f.has(Integral):
                return None
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
        from sympy.solvers.inequalities import _solve_inequality
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


def _complete_the_square_in_denom(f, s):
    from sympy.simplify.radsimp import fraction
    [n, d] = fraction(f)
    if d.is_polynomial(s):
        cf = d.as_poly(s).all_coeffs()
        if len(cf) == 3:
            a, b, c = cf
            d = a*((s+b/(2*a))**2+c/a-(b/(2*a))**2)
    return n/d


@cacheit
def _inverse_laplace_build_rules():
    """
    This is an internal helper function that returns the table of inverse
    Laplace transform rules in terms of the time variable `t` and the
    frequency variable `s`.  It is used by `_inverse_laplace_apply_rules`.
    """
    s = Dummy('s')
    t = Dummy('t')
    a = Wild('a', exclude=[s])
    b = Wild('b', exclude=[s])
    c = Wild('c', exclude=[s])

    debug('_inverse_laplace_build_rules is building rules')

    def _frac(f, s):
        try:
            return f.factor(s)
        except PolynomialError:
            return f

    def same(f): return f
    # This list is sorted according to the prep function needed.
    _ILT_rules = [
        (a/s, a, S.true, same, 1),
        (b*(s+a)**(-c), t**(c-1)*exp(-a*t)/gamma(c), c > 0, same, 1),
        (1/(s**2+a**2)**2, (sin(a*t) - a*t*cos(a*t))/(2*a**3),
         S.true, same, 1),
        # The next two rules must be there in that order. For the second
        # one, the condition would be a != 0 or, respectively, to take the
        # limit a -> 0 after the transform if a == 0. It is much simpler if
        # the case a == 0 has its own rule.
        (1/(s**b), t**(b - 1)/gamma(b), S.true, same, 1),
        (1/(s*(s+a)**b), lowergamma(b, a*t)/(a**b*gamma(b)),
         S.true, same, 1)
    ]
    return _ILT_rules, s, t


def _inverse_laplace_apply_simple_rules(f, s, t):
    """
    Helper function for the class InverseLaplaceTransform.
    """
    if f == 1:
        debug('_inverse_laplace_apply_simple_rules match:')
        debugf('      f:    %s', (1,))
        debugf('      rule: 1 o---o DiracDelta(%s)', (t,))
        return DiracDelta(t), S.true

    _ILT_rules, s_, t_ = _inverse_laplace_build_rules()
    _prep = ''
    fsubs = f.subs({s: s_})

    for s_dom, t_dom, check, prep, fac in _ILT_rules:
        if _prep != (prep, fac):
            _F = prep(fsubs*fac)
            _prep = (prep, fac)
        ma = _F.match(s_dom)
        if ma:
            try:
                c = check.xreplace(ma)
            except TypeError:
                continue
            if c == S.true:
                debug('_inverse_laplace_apply_simple_rules match:')
                debugf('      f:    %s', (f,))
                debugf('      rule: %s o---o %s', (s_dom, t_dom))
                debugf('      ma:   %s', (ma,))
                return Heaviside(t)*t_dom.xreplace(ma).subs({t_: t}), S.true

    return None


def _inverse_laplace_time_shift(F, s, t, plane):
    """
    Helper function for the class InverseLaplaceTransform.
    """
    a = Wild('a', exclude=[s])
    g = Wild('g')

    if not F.has(s):
        return F*DiracDelta(t), S.true
    ma1 = F.match(exp(a*s))
    if ma1:
        if ma1[a].is_negative:
            debug('_inverse_laplace_time_shift match:')
            debugf('      f:    %s', (F,))
            debug('      rule: exp(-a*s) o---o DiracDelta(t-a)')
            debugf('      ma:   %s', (ma1,))
            return DiracDelta(t+ma1[a]), S.true
        else:
            debug('_inverse_laplace_time_shift match: negative time shift')
            return InverseLaplaceTransform(F, s, t, plane), S.true

    ma1 = F.match(exp(a*s)*g)
    if ma1:
        if ma1[a].is_negative:
            debug('_inverse_laplace_time_shift match:')
            debugf('      f:    %s', (F,))
            debug('      rule: exp(-a*s)*F(s) o---o Heaviside(t-a)*f(t-a)')
            debugf('      ma:   %s', (ma1,))
            return _inverse_laplace_transform(ma1[g], s, t+ma1[a], plane)
        else:
            debug('_inverse_laplace_time_shift match: negative time shift')
            return InverseLaplaceTransform(F, s, t, plane), S.true
    return None


def _inverse_laplace_time_diff(F, s, t, plane):
    """
    Helper function for the class InverseLaplaceTransform.
    """
    n = Wild('n', exclude=[s])
    g = Wild('g')

    ma1 = F.match(s**n*g)
    if ma1 and ma1[n].is_integer and ma1[n].is_positive:
        debug('_inverse_laplace_time_diff match:')
        debugf('      f:    %s', (F,))
        debug('      rule: s**n*F(s) o---o diff(f(t), t, n)')
        debugf('      ma:   %s', (ma1,))
        r, c = _inverse_laplace_transform(ma1[g], s, t, plane)
        r = r.replace(Heaviside(t), 1)
        if r.has(InverseLaplaceTransform):
            return diff(r, t, ma1[n]), c
        else:
            return Heaviside(t)*diff(r, t, ma1[n]), c
    return None


def _inverse_laplace_apply_prog_rules(F, s, t, plane):
    """
    Helper function for the class InverseLaplaceTransform.
    """
    prog_rules = [_inverse_laplace_time_shift,
                  _inverse_laplace_time_diff]

    for p_rule in prog_rules:
        if (r := p_rule(F, s, t, plane)) is not None:
            return r
    return None


def _inverse_laplace_expand(fn, s, t, plane):
    """
    Helper function for the class InverseLaplaceTransform.
    """
    if fn.is_Add:
        return None
    r = expand(fn, deep=False)
    if r.is_Add:
        return _inverse_laplace_transform(r, s, t, plane)
    r = expand_mul(fn)
    if r.is_Add:
        return _inverse_laplace_transform(r, s, t, plane)
    r = expand(fn)
    if r.is_Add:
        return _inverse_laplace_transform(r, s, t, plane)
    if fn.is_rational_function(s):
        r = fn.apart(s).doit()
    if r.is_Add:
        return _inverse_laplace_transform(r, s, t, plane)
    return None


def _inverse_laplace_rational(fn, s, t, plane, simplify):
    """
    Helper function for the class InverseLaplaceTransform.
    """
    debugf('[ILT _i_l_r] (%s, %s, %s)', (fn, s, t))
    x_ = symbols('x_')
    f = fn.apart(s)
    terms = Add.make_args(f)
    terms_t = []
    conditions = [S.true]
    for term in terms:
        [n, d] = term.as_numer_denom()
        dc = d.as_poly(s).all_coeffs()
        dc_lead = dc[0]
        dc = [x/dc_lead for x in dc]
        nc = [x/dc_lead for x in n.as_poly(s).all_coeffs()]
        if len(dc) == 1:
            r = nc[0]*DiracDelta(t)
            terms_t.append(r)
        elif len(dc) == 2:
            r = nc[0]*exp(-dc[1]*t)
            terms_t.append(Heaviside(t)*r)
        elif len(dc) == 3:
            a = dc[1]/2
            b = (dc[2]-a**2).factor()
            if len(nc) == 1:
                nc = [S.Zero] + nc
            l, m = tuple(nc)
            if b == 0:
                r = (m*t+l*(1-a*t))*exp(-a*t)
            else:
                hyp = False
                if b.is_negative:
                    b = -b
                    hyp = True
                b2 = list(roots(x_**2-b, x_).keys())[0]
                bs = sqrt(b).simplify()
                if hyp:
                    r = (
                        l*exp(-a*t)*cosh(b2*t) + (m-a*l) /
                        bs*exp(-a*t)*sinh(bs*t))
                else:
                    r = l*exp(-a*t)*cos(b2*t) + (m-a*l)/bs*exp(-a*t)*sin(bs*t)
            terms_t.append(Heaviside(t)*r)
        else:
            ft, cond = _inverse_laplace_transform(
                fn, s, t, plane, simplify=True, dorational=False)
            terms_t.append(ft)
            conditions.append(cond)

    result = Add(*terms_t)
    if simplify:
        result = result.simplify(doit=False)
    debugf('[ILT _i_l_r]   returns %s', (result,))
    return result, And(*conditions)


def _inverse_laplace_transform(
        fn, s_, t_, plane, simplify=True, dorational=True):
    """
    Front-end function of the inverse Laplace transform. It tries to apply all
    known rules recursively.  If everything else fails, it tries to integrate.
    """
    terms = Add.make_args(fn)
    terms_t = []
    conditions = []

    debugf('[ILT _i_l_t] (%s, %s, %s)', (fn, s_, t_))

    for term in terms:
        k, f = term.as_independent(s_, as_Add=False)
        if (
                dorational and term.is_rational_function(s_) and
                (
                    r := _inverse_laplace_rational(
                        f, s_, t_, plane, simplify)) is not None):
            pass
        elif (r := _inverse_laplace_apply_simple_rules(f, s_, t_)) is not None:
            pass
        elif (r := _inverse_laplace_expand(f, s_, t_, plane)) is not None:
            pass
        elif (
                (r := _inverse_laplace_apply_prog_rules(f, s_, t_, plane))
                is not None):
            pass
        elif any(undef.has(s_) for undef in f.atoms(AppliedUndef)):
            # If there are undefined functions f(t) then integration is
            # unlikely to do anything useful so we skip it and given an
            # unevaluated LaplaceTransform.
            r = (InverseLaplaceTransform(f, s_, t_, plane), S.true)
        elif (
                r := _inverse_laplace_transform_integration(
                    f, s_, t_, plane, simplify=simplify)) is not None:
            pass
        else:
            r = (InverseLaplaceTransform(f, s_, t_, plane), S.true)
        (ri_, ci_) = r
        terms_t.append(k*ri_)
        conditions.append(ci_)

    result = Add(*terms_t)
    if simplify:
        result = result.simplify(doit=False)
    condition = And(*conditions)

    return result, condition


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
        return _inverse_laplace_transform_integration(
            F, s, t, self.fundamental_plane, **hints)

    def _as_integral(self, F, s, t):
        c = self.__class__._c
        return (
            Integral(exp(s*t)*F, (s, c - S.ImaginaryUnit*S.Infinity,
                                  c + S.ImaginaryUnit*S.Infinity)) /
            (2*S.Pi*S.ImaginaryUnit))

    def doit(self, **hints):
        """
        Try to evaluate the transform in closed form.

        Explanation
        ===========

        Standard hints are the following:
        - ``noconds``:  if True, do not return convergence conditions. The
        default setting is `True`.
        - ``simplify``: if True, it simplifies the final result. The
        default setting is `False`.
        """
        _noconds = hints.get('noconds', True)
        _simplify = hints.get('simplify', False)

        debugf('[ILT doit] (%s, %s, %s)', (self.function,
                                           self.function_variable,
                                           self.transform_variable))

        s_ = self.function_variable
        t_ = self.transform_variable
        fn = self.function
        plane = self.fundamental_plane

        r = _inverse_laplace_transform(fn, s_, t_, plane, simplify=_simplify)

        if _noconds:
            return r[0]
        else:
            return r


def inverse_laplace_transform(F, s, t, plane=None, **hints):
    r"""
    Compute the inverse Laplace transform of `F(s)`, defined as

    .. math ::
        f(t) = \frac{1}{2\pi i} \int_{c-i\infty}^{c+i\infty} e^{st}
        F(s) \mathrm{d}s,

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

    laplace_transform
    hankel_transform, inverse_hankel_transform
    """
    if isinstance(F, MatrixBase) and hasattr(F, 'applyfunc'):
        return F.applyfunc(
            lambda Fij: inverse_laplace_transform(Fij, s, t, plane, **hints))
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
