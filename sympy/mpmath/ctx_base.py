from operator import gt, lt

from functions.functions import SpecialFunctions
from functions.rszeta import RSCache
from calculus.quadrature import QuadratureMethods
from calculus.calculus import CalculusMethods
from calculus.optimization import OptimizationMethods
from calculus.odes import ODEMethods
from matrices.matrices import MatrixMethods
from matrices.calculus import MatrixCalculusMethods
from matrices.linalg import LinearAlgebraMethods
from identification import IdentificationMethods
from visualization import VisualizationMethods

import libmp

class Context(object):
    pass

class StandardBaseContext(Context,
    SpecialFunctions,
    RSCache,
    QuadratureMethods,
    CalculusMethods,
    MatrixMethods,
    MatrixCalculusMethods,
    LinearAlgebraMethods,
    IdentificationMethods,
    OptimizationMethods,
    ODEMethods,
    VisualizationMethods):

    NoConvergence = libmp.NoConvergence
    ComplexResult = libmp.ComplexResult

    def __init__(ctx):
        ctx._aliases = {}
        # Call those that need preinitialization (e.g. for wrappers)
        SpecialFunctions.__init__(ctx)
        RSCache.__init__(ctx)
        QuadratureMethods.__init__(ctx)
        CalculusMethods.__init__(ctx)
        MatrixMethods.__init__(ctx)

    def _init_aliases(ctx):
        for alias, value in ctx._aliases.items():
            try:
                setattr(ctx, alias, getattr(ctx, value))
            except AttributeError:
                pass

    _fixed_precision = False

    # XXX
    verbose = False

    def warn(ctx, msg):
        print "Warning:", msg

    def bad_domain(ctx, msg):
        raise ValueError(msg)

    def _re(ctx, x):
        if hasattr(x, "real"):
            return x.real
        return x

    def _im(ctx, x):
        if hasattr(x, "imag"):
            return x.imag
        return ctx.zero

    def chop(ctx, x, tol=None):
        """
        Chops off small real or imaginary parts, or converts
        numbers close to zero to exact zeros. The input can be a
        single number or an iterable::

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> chop(5+1e-10j, tol=1e-9)
            mpf('5.0')
            >>> nprint(chop([1.0, 1e-20, 3+1e-18j, -4, 2]))
            [1.0, 0.0, 3.0, -4.0, 2.0]

        The tolerance defaults to ``100*eps``.
        """
        if tol is None:
            tol = 100*ctx.eps
        try:
            x = ctx.convert(x)
            absx = abs(x)
            if abs(x) < tol:
                return ctx.zero
            if ctx._is_complex_type(x):
                if abs(x.imag) < min(tol, absx*tol):
                    return x.real
                if abs(x.real) < min(tol, absx*tol):
                    return ctx.mpc(0, x.imag)
        except TypeError:
            if isinstance(x, ctx.matrix):
                return x.apply(lambda a: ctx.chop(a, tol))
            if hasattr(x, "__iter__"):
                return [ctx.chop(a, tol) for a in x]
        return x

    def almosteq(ctx, s, t, rel_eps=None, abs_eps=None):
        r"""
        Determine whether the difference between `s` and `t` is smaller
        than a given epsilon, either relatively or absolutely.

        Both a maximum relative difference and a maximum difference
        ('epsilons') may be specified. The absolute difference is
        defined as `|s-t|` and the relative difference is defined
        as `|s-t|/\max(|s|, |t|)`.

        If only one epsilon is given, both are set to the same value.
        If none is given, both epsilons are set to `2^{-p+m}` where
        `p` is the current working precision and `m` is a small
        integer. The default setting typically allows :func:`almosteq`
        to be used to check for mathematical equality
        in the presence of small rounding errors.

        **Examples**

            >>> from mpmath import *
            >>> mp.dps = 15
            >>> almosteq(3.141592653589793, 3.141592653589790)
            True
            >>> almosteq(3.141592653589793, 3.141592653589700)
            False
            >>> almosteq(3.141592653589793, 3.141592653589700, 1e-10)
            True
            >>> almosteq(1e-20, 2e-20)
            True
            >>> almosteq(1e-20, 2e-20, rel_eps=0, abs_eps=0)
            False

        """
        t = ctx.convert(t)
        if abs_eps is None and rel_eps is None:
            rel_eps = abs_eps = ctx.ldexp(1, -ctx.prec+4)
        if abs_eps is None:
            abs_eps = rel_eps
        elif rel_eps is None:
            rel_eps = abs_eps
        diff = abs(s-t)
        if diff <= abs_eps:
            return True
        abss = abs(s)
        abst = abs(t)
        if abss < abst:
            err = diff/abst
        else:
            err = diff/abss
        return err <= rel_eps

    def arange(ctx, *args):
        r"""
        This is a generalized version of Python's :func:`range` function
        that accepts fractional endpoints and step sizes and
        returns a list of ``mpf`` instances. Like :func:`range`,
        :func:`arange` can be called with 1, 2 or 3 arguments:

        ``arange(b)``
            `[0, 1, 2, \ldots, x]`
        ``arange(a, b)``
            `[a, a+1, a+2, \ldots, x]`
        ``arange(a, b, h)``
            `[a, a+h, a+h, \ldots, x]`

        where `b-1 \le x < b` (in the third case, `b-h \le x < b`).

        Like Python's :func:`range`, the endpoint is not included. To
        produce ranges where the endpoint is included, :func:`linspace`
        is more convenient.

        **Examples**

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> arange(4)
            [mpf('0.0'), mpf('1.0'), mpf('2.0'), mpf('3.0')]
            >>> arange(1, 2, 0.25)
            [mpf('1.0'), mpf('1.25'), mpf('1.5'), mpf('1.75')]
            >>> arange(1, -1, -0.75)
            [mpf('1.0'), mpf('0.25'), mpf('-0.5')]

        """
        if not len(args) <= 3:
            raise TypeError('arange expected at most 3 arguments, got %i'
                            % len(args))
        if not len(args) >= 1:
            raise TypeError('arange expected at least 1 argument, got %i'
                            % len(args))
        # set default
        a = 0
        dt = 1
        # interpret arguments
        if len(args) == 1:
            b = args[0]
        elif len(args) >= 2:
            a = args[0]
            b = args[1]
        if len(args) == 3:
            dt = args[2]
        a, b, dt = ctx.mpf(a), ctx.mpf(b), ctx.mpf(dt)
        assert a + dt != a, 'dt is too small and would cause an infinite loop'
        # adapt code for sign of dt
        if a > b:
            if dt > 0:
                return []
            op = gt
        else:
            if dt < 0:
                return []
            op = lt
        # create list
        result = []
        i = 0
        t = a
        while 1:
            t = a + dt*i
            i += 1
            if op(t, b):
                result.append(t)
            else:
                break
        return result

    def linspace(ctx, *args, **kwargs):
        """
        ``linspace(a, b, n)`` returns a list of `n` evenly spaced
        samples from `a` to `b`. The syntax ``linspace(mpi(a,b), n)``
        is also valid.

        This function is often more convenient than :func:`arange`
        for partitioning an interval into subintervals, since
        the endpoint is included::

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> linspace(1, 4, 4)
            [mpf('1.0'), mpf('2.0'), mpf('3.0'), mpf('4.0')]
            >>> linspace(mpi(1,4), 4)
            [mpf('1.0'), mpf('2.0'), mpf('3.0'), mpf('4.0')]

        You may also provide the keyword argument ``endpoint=False``::

            >>> linspace(1, 4, 4, endpoint=False)
            [mpf('1.0'), mpf('1.75'), mpf('2.5'), mpf('3.25')]

        """
        if len(args) == 3:
            a = ctx.mpf(args[0])
            b = ctx.mpf(args[1])
            n = int(args[2])
        elif len(args) == 2:
            assert hasattr(args[0], '_mpi_')
            a = args[0].a
            b = args[0].b
            n = int(args[1])
        else:
            raise TypeError('linspace expected 2 or 3 arguments, got %i' \
                            % len(args))
        if n < 1:
            raise ValueError('n must be greater than 0')
        if not 'endpoint' in kwargs or kwargs['endpoint']:
            if n == 1:
                return [ctx.mpf(a)]
            step = (b - a) / ctx.mpf(n - 1)
            y = [i*step + a for i in xrange(n)]
            y[-1] = b
        else:
            step = (b - a) / ctx.mpf(n)
            y = [i*step + a for i in xrange(n)]
        return y

    def cos_sin(ctx, z, **kwargs):
        return ctx.cos(z, **kwargs), ctx.sin(z, **kwargs)

    def _default_hyper_maxprec(ctx, p):
        return int(1000 * p**0.25 + 4*p)

    _gcd = staticmethod(libmp.gcd)
    list_primes = staticmethod(libmp.list_primes)
    bernfrac = staticmethod(libmp.bernfrac)
    moebius = staticmethod(libmp.moebius)
    _ifac = staticmethod(libmp.ifac)
    _eulernum = staticmethod(libmp.eulernum)

    def sum_accurately(ctx, terms, check_step=1):
        prec = ctx.prec
        try:
            extraprec = 10
            while 1:
                ctx.prec = prec + extraprec + 5
                max_mag = ctx.ninf
                s = ctx.zero
                k = 0
                for term in terms():
                    s += term
                    if (not k % check_step) and term:
                        term_mag = ctx.mag(term)
                        max_mag = max(max_mag, term_mag)
                        sum_mag = ctx.mag(s)
                        if sum_mag - term_mag > ctx.prec:
                            break
                    k += 1
                cancellation = max_mag - sum_mag
                if cancellation != cancellation:
                    break
                if cancellation < extraprec or ctx._fixed_precision:
                    break
                extraprec += min(ctx.prec, cancellation)
            return s
        finally:
            ctx.prec = prec

    def power(ctx, x, y):
        return ctx.convert(x) ** ctx.convert(y)

    def _zeta_int(ctx, n):
        return ctx.zeta(n)
