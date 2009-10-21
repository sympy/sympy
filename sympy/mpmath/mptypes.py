"""
This module defines the mpf, mpc classes, and standard functions for
operating with them.
"""
__docformat__ = 'plaintext'

import re

from settings import (MP_BASE, MP_ZERO, MP_ONE, int_types, repr_dps,
    round_floor, round_ceiling, dps_to_prec, round_nearest, prec_to_dps)

from libmpf import (
    ComplexResult, to_pickable, from_pickable, normalize,
    from_int, from_float, from_str, to_int, to_float, to_str,
    from_rational, from_man_exp,
    fone, fzero, finf, fninf, fnan,
    mpf_abs, mpf_pos, mpf_neg, mpf_add, mpf_sub, mpf_mul, mpf_mul_int,
    mpf_div, mpf_rdiv_int, mpf_pow_int, mpf_mod,
    mpf_eq, mpf_cmp, mpf_lt, mpf_gt, mpf_le, mpf_ge,
    mpf_hash, mpf_rand,
    mpf_sum
)

from libmpc import (
    mpc_to_str,
    mpc_to_complex, mpc_hash, mpc_pos, mpc_is_nonzero, mpc_neg, mpc_conjugate,
    mpc_abs, mpc_add, mpc_add_mpf, mpc_sub, mpc_sub_mpf, mpc_mul, mpc_mul_mpf,
    mpc_mul_int, mpc_div, mpc_div_mpf, mpc_pow, mpc_pow_mpf, mpc_pow_int,
    mpc_mpf_div
)

from libelefun import mpf_pow

from libmpi import (
    mpi_mid, mpi_delta, mpi_str,
    mpi_abs, mpi_pos, mpi_neg, mpi_add, mpi_sub,
    mpi_mul, mpi_div, mpi_pow_int, mpi_pow
)

from string import strip
from operator import gt, lt


new = object.__new__

class PrecisionManager:

    def __init__(self, mp, precfun, dpsfun, normalize_output=False):
        self.mp = mp
        self.precfun = precfun
        self.dpsfun = dpsfun
        self.normalize_output = normalize_output

    def __call__(self, f):
        def g(*args, **kwargs):
            orig = self.mp.prec
            try:
                if self.precfun:
                    self.mp.prec = self.precfun(self.mp.prec)
                else:
                    self.mp.dps = self.dpsfun(self.mp.dps)
                if self.normalize_output:
                    v = f(*args, **kwargs)
                    if type(v) is tuple:
                        return tuple([+a for a in v])
                    return +v
                else:
                    return f(*args, **kwargs)
            finally:
                self.mp.prec = orig
        g.__name__ = f.__name__
        g.__doc__ = f.__doc__
        return g

    def __enter__(self):
        self.origp = self.mp.prec
        if self.precfun:
            self.mp.prec = self.precfun(self.mp.prec)
        else:
            self.mp.dps = self.dpsfun(self.mp.dps)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.mp.prec = self.origp
        return False

class Context(object):
    pass


class MultiPrecisionArithmetic(Context):
    """
    Context for multiprecision arithmetic with a global precision.
    """

    def __init__(ctx):
        # Settings
        ctx._prec_rounding = [53, round_nearest]
        ctx.trap_complex = False
        ctx.mpf = type('mpf', (_mpf,), {})
        ctx.mpc = type('mpc', (_mpc,), {})
        ctx.mpi = type('mpi', (_mpi,), {})
        ctx.constant = type('constant', (_constant,), {})
        ctx.types = [ctx.mpf, ctx.mpc, ctx.mpi, ctx.constant]
        # For fast access
        ctx.mpf._ctxdata = [ctx.mpf, new, ctx._prec_rounding]
        ctx.mpc._ctxdata = [ctx.mpc, new, ctx._prec_rounding]
        ctx.mpi._ctxdata = [ctx.mpi, new, ctx._prec_rounding]
        ctx.constant._ctxdata = [ctx.mpf, new, ctx._prec_rounding]
        ctx.constant.context = ctx
        ctx.mpf.context = ctx
        ctx.mpc.context = ctx
        ctx.mpi.context = ctx
        # Predefined data
        ctx._create_constants({})

    # pi, etc
    _constants = []

    def _create_constants(ctx, namespace):
        ctx.one = ctx.make_mpf(fone)
        ctx.zero = ctx.make_mpf(fzero)
        ctx.inf = ctx.make_mpf(finf)
        ctx.ninf = ctx.make_mpf(fninf)
        ctx.nan = ctx.make_mpf(fnan)
        ctx.j = ctx.make_mpc((fzero,fone))
        eps = ctx.constant(lambda prec, rnd: (0, MP_ONE, 1-prec, 1),
            "epsilon of working precision")
        ctx.eps = eps
        import function_docs
        for name, func, descr in ctx._constants:
            doc = function_docs.__dict__.get(name, descr)
            const_cls = type("_" + name, (ctx.constant,), {'__doc__':doc})
            const_inst = const_cls(func, descr)
            setattr(ctx, name, const_inst)
            namespace[name] = const_inst

    def clone(ctx):
        """
        Create a copy of the context, with the same working precision.
        """
        a = MultiPrecisionArithmetic()
        a.prec = ctx.prec
        return a

    # Several helper methods
    # TODO: add more of these, make consistent, write docstrings, ...

    def is_real_type(ctx, x):
        return hasattr(x, '_mpf_')

    def is_complex_type(ctx, x):
        return hasattr(x, '_mpc_')

    def make_mpf(ctx, v):
        a = new(ctx.mpf)
        a._mpf_ = v
        return a

    def make_mpc(ctx, v):
        a = new(ctx.mpc)
        a._mpc_ = v
        return a

    def make_mpi(ctx, v):
        a = new(ctx.mpi)
        a._mpi_ = v
        return a

    def isnpint(ctx, x):
        if not x:
            return True
        if hasattr(x, '_mpf_'):
            sign, man, exp, bc = x._mpf_
            return sign and exp >= 0
        if hasattr(x, '_mpc_'):
            return not x.imag and ctx.isnpint(x.real)

    def bad_domain(ctx, msg):
        raise ValueError(msg)

    def __str__(ctx):
        lines = ["Mpmath settings:",
            ("  mp.prec = %s" % ctx.prec).ljust(30) + "[default: 53]",
            ("  mp.dps = %s" % ctx.dps).ljust(30) + "[default: 15]",
            ("  mp.trap_complex = %s" % ctx.trap_complex).ljust(30) + "[default: False]",
        ]
        return "\n".join(lines)

    def default(ctx):
        ctx._prec = ctx._prec_rounding[0] = 53
        ctx._dps = 15
        ctx.trap_complex = False

    def set_prec(ctx, n):
        ctx._prec = ctx._prec_rounding[0] = max(1, int(n))
        ctx._dps = prec_to_dps(n)

    def set_dps(ctx, n):
        ctx._prec = ctx._prec_rounding[0] = dps_to_prec(n)
        ctx._dps = max(1, int(n))

    prec = property(lambda ctx: ctx._prec, set_prec)
    dps = property(lambda ctx: ctx._dps, set_dps)

    @property
    def repr_digits(ctx):
        return repr_dps(ctx._prec)

    @property
    def str_digits(ctx):
        return ctx._dps

    verbose = False

    def extraprec(ctx, n, normalize_output=False):
        """
        The block

            with extraprec(n):
                <code>

        increases the precision n bits, executes <code>, and then
        restores the precision.

        extraprec(n)(f) returns a decorated version of the function f
        that increases the working precision by n bits before execution,
        and restores the parent precision afterwards. With
        normalize_output=True, it rounds the return value to the parent
        precision.
        """
        return PrecisionManager(ctx, lambda p: p + n, None, normalize_output)

    def extradps(ctx, n, normalize_output=False):
        """
        This function is analogous to extraprec (see documentation)
        but changes the decimal precision instead of the number of bits.
        """
        return PrecisionManager(ctx, None, lambda d: d + n, normalize_output)

    def workprec(ctx, n, normalize_output=False):
        """
        The block

            with workprec(n):
                <code>

        sets the precision to n bits, executes <code>, and then restores
        the precision.

        workprec(n)(f) returns a decorated version of the function f
        that sets the precision to n bits before execution,
        and restores the precision afterwards. With normalize_output=True,
        it rounds the return value to the parent precision.
        """
        return PrecisionManager(ctx, lambda p: n, None, normalize_output)

    def workdps(ctx, n, normalize_output=False):
        """
        This function is analogous to workprec (see documentation)
        but changes the decimal precision instead of the number of bits.
        """
        return PrecisionManager(ctx, None, lambda d: n, normalize_output)

    def nstr(ctx, x, n=6, **kwargs):
        """
        Convert an ``mpf``, ``mpc`` or ``mpi`` to a decimal string literal with *n*
        significant digits. The small default value for *n* is chosen to
        make this function useful for printing collections of numbers
        (lists, matrices, etc).

        If *x* is an ``mpi``, there are some extra options, notably *mode*, which
        can be 'brackets', 'diff', 'plusminus' or 'percent'. See ``mpi_to_str`` for
        a more complete documentation.

        If *x* is a list or tuple, :func:`nstr` is applied recursively
        to each element. For unrecognized classes, :func:`nstr`
        simply returns ``str(x)``.

        The companion function :func:`nprint` prints the result
        instead of returning it.

            >>> from mpmath import *
            >>> nstr([+pi, ldexp(1,-500)])
            '[3.14159, 3.05494e-151]'
            >>> nprint([+pi, ldexp(1,-500)])
            [3.14159, 3.05494e-151]

        """
        if isinstance(x, list):
            return "[%s]" % (", ".join(nstr(c, n) for c in x))
        if isinstance(x, tuple):
            return "(%s)" % (", ".join(nstr(c, n) for c in x))
        if hasattr(x, '_mpf_'):
            return to_str(x._mpf_, n)
        if hasattr(x, '_mpc_'):
            return "(" + mpc_to_str(x._mpc_, n)  + ")"
        if isinstance(x, basestring):
            return repr(x)
        from matrices import matrix
        if isinstance(x, matrix):
            return x.__nstr__(n)
        if hasattr(x, '_mpi_'):
            return ctx.mpi_to_str(x, n, **kwargs)
        return str(x)

    def nprint(ctx, x, n=6, **kwargs):
        """
        Equivalent to ``print nstr(x, n)``.
        """
        print ctx.nstr(x, n, **kwargs)

    def convert(ctx, x, strings=True):
        """
        Converts *x* to an ``mpf``, ``mpc`` or ``mpi``. If *x* is of type ``mpf``,
        ``mpc``, ``int``, ``float``, ``complex``, the conversion
        will be performed losslessly.

        If *x* is a string, the result will be rounded to the present
        working precision. Strings representing fractions or complex
        numbers are permitted.

            >>> from mpmath import *
            >>> mp.dps = 15
            >>> mpmathify(3.5)
            mpf('3.5')
            >>> mpmathify('2.1')
            mpf('2.1000000000000001')
            >>> mpmathify('3/4')
            mpf('0.75')
            >>> mpmathify('2+3j')
            mpc(real='2.0', imag='3.0')

        """
        if type(x) in ctx.types: return x
        if isinstance(x, int_types): return ctx.make_mpf(from_int(x))
        if isinstance(x, float): return ctx.make_mpf(from_float(x))
        if isinstance(x, complex): return ctx.mpc(x)
        prec, rounding = ctx._prec_rounding
        if strings and isinstance(x, basestring):
            try:
                _mpf_ = from_str(x, prec, rounding)
                return ctx.make_mpf(_mpf_)
            except Exception, e:
                if '/' in x:
                    fract = x.split('/')
                    assert len(fract) == 2
                    return ctx.convert(fract[0]) / ctx.convert(fract[1])
                if 'j' in x.lower():
                    x = x.lower().replace(' ', '')
                    match = get_complex.match(x)
                    re = match.group('re')
                    if not re:
                        re = 0
                    im = match.group('im').rstrip('j')
                    return ctx.mpc(ctx.convert(re), ctx.convert(im))
                if '[' in x or '(' in x or '+-' in x:
                    # XXX
                    return ctx.mpi_from_str(x)
                raise e
        if hasattr(x, '_mpf_'): return ctx.make_mpf(x._mpf_)
        if hasattr(x, '_mpc_'): return ctx.make_mpc(x._mpc_)
        if hasattr(x, '_mpmath_'):
            return ctx.convert(x._mpmath_(*prec_rounding))
        raise TypeError("cannot create mpf from " + repr(x))

    def isnan(ctx, x):
        """
        For an ``mpf`` *x*, determines whether *x* is not-a-number (nan)::

            >>> from mpmath import *
            >>> isnan(nan), isnan(3)
            (True, False)

        """
        if not hasattr(x, '_mpf_'):
            return False
        return x._mpf_ == fnan

    def isinf(ctx, x):
        """
        For an ``mpf`` *x*, determines whether *x* is infinite::

            >>> from mpmath import *
            >>> isinf(inf), isinf(-inf), isinf(3)
            (True, True, False)

        """
        if not hasattr(x, '_mpf_'):
            return False
        return x._mpf_ in (finf, fninf)

    def isint(ctx, x):
        """
        For an ``mpf`` *x*, or any type that can be converted
        to ``mpf``, determines whether *x* is exactly
        integer-valued::

            >>> from mpmath import *
            >>> isint(3), isint(mpf(3)), isint(3.2)
            (True, True, False)

        """
        if isinstance(x, int_types):
            return True
        try:
            x = ctx.convert(x)
        except:
            return False
        if hasattr(x, '_mpf_'):
            if ctx.isnan(x) or ctx.isinf(x):
                return False
            return x == int(x)
        return False

    def chop(ctx, x, tol=None):
        """
        Chops off small real or imaginary parts, or converts
        numbers close to zero to exact zeros. The input can be a
        single number or an iterable::

            >>> from mpmath import *
            >>> mp.dps = 15
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
            if ctx.is_complex_type(x):
                if abs(x.imag) < min(tol, absx*tol):
                    return x.real
                if abs(x.real) < min(tol, absx*tol):
                    return ctx.mpc(0, x.imag)
        except TypeError:
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
            rel_eps = abs_eps = ctx.make_mpf((0, 1, -ctx.prec+4, 1))
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

    def fsum(ctx, terms, absolute=False, squared=False):
        """
        Calculates a sum containing a finite number of terms (for infinite
        series, see :func:`nsum`). The terms will be converted to
        mpmath numbers. For len(terms) > 2, this function is generally
        faster and produces more accurate results than the builtin
        Python function :func:`sum`.

            >>> from mpmath import *
            >>> mp.dps = 15
            >>> fsum([1, 2, 0.5, 7])
            mpf('10.5')

        With squared=True each term is squared, and with absolute=True
        the absolute value of each term is used.
        """
        prec, rnd = ctx._prec_rounding
        real = []
        imag = []
        other = 0
        for term in terms:
            reval = imval = 0
            if hasattr(term, "_mpf_"):
                reval = term._mpf_
            elif hasattr(term, "_mpc_"):
                reval, imval = term._mpc_
            else:
                term = ctx.convert(term)
                if hasattr(term, "_mpf_"):
                    reval = term._mpf_
                elif hasattr(term, "_mpc_"):
                    reval, imval = term._mpc_
                else:
                    if absolute: term = ctx.absmax(term)
                    if squared: term = term**2
                    other += term
                    continue
            if imval:
                if squared:
                    if absolute:
                        real.append(mpf_mul(reval,reval))
                        real.append(mpf_mul(imval,imval))
                    else:
                        reval, imval = mpc_pow_int((reval,imval),2,prec+10)
                        real.append(reval)
                        imag.append(imval)
                elif absolute:
                    real.append(mpc_abs((reval,imval), prec))
                else:
                    real.append(reval)
                    imag.append(imval)
            else:
                if squared:
                    reval = mpf_mul(reval, reval)
                elif absolute:
                    reval = mpf_abs(reval)
                real.append(reval)
        s = mpf_sum(real, prec, rnd, absolute)
        if imag:
            s = ctx.make_mpc((s, mpf_sum(imag, prec, rnd)))
        else:
            s = ctx.make_mpf(s)
        if other is 0:
            return s
        else:
            return s + other

    def fdot(ctx, A, B=None):
        r"""
        Computes the dot product of the iterables `A` and `B`,

        .. math ::

            \sum_{k=0} A_k B_k.

        Alternatively, :func:`fdot` accepts a single iterable of pairs.
        In other words, ``fdot(A,B)`` and ``fdot(zip(A,B))`` are equivalent.

        The elements are automatically converted to mpmath numbers.

        Examples::

            >>> from mpmath import *
            >>> mp.dps = 15
            >>> A = [2, 1.5, 3]
            >>> B = [1, -1, 2]
            >>> fdot(A, B)
            mpf('6.5')
            >>> zip(A, B)
            [(2, 1), (1.5, -1), (3, 2)]
            >>> fdot(_)
            mpf('6.5')

        """
        if B:
            A = zip(A, B)
        prec, rnd = ctx._prec_rounding
        real = []
        imag = []
        other = 0
        hasattr_ = hasattr
        types = (ctx.mpf, ctx.mpc)
        for a, b in A:
            if type(a) not in types: a = ctx.convert(a)
            if type(b) not in types: b = ctx.convert(b)
            a_real = hasattr_(a, "_mpf_")
            b_real = hasattr_(b, "_mpf_")
            if a_real and b_real:
                real.append(mpf_mul(a._mpf_, b._mpf_))
                continue
            a_complex = hasattr_(a, "_mpc_")
            b_complex = hasattr_(b, "_mpc_")
            if a_real and b_complex:
                aval = a._mpf_
                bre, bim = b._mpc_
                real.append(mpf_mul(aval, bre))
                imag.append(mpf_mul(aval, bim))
            elif b_real and a_complex:
                are, aim = a._mpc_
                bval = b._mpf_
                real.append(mpf_mul(are, bval))
                imag.append(mpf_mul(aim, bval))
            elif a_complex and b_complex:
                re, im = mpc_mul(a._mpc_, b._mpc_, prec+20)
                real.append(re)
                imag.append(im)
            else:
                other += a*b
        s = mpf_sum(real, prec, rnd)
        if imag:
            s = ctx.make_mpc((s, mpf_sum(imag, prec, rnd)))
        else:
            s = ctx.make_mpf(s)
        if other is 0:
            return s
        else:
            return s + other

    def fprod(ctx, factors):
        r"""
        Calculates a product containing a finite number of factors (for
        infinite products, see :func:`nprod`). The factors will be
        converted to mpmath numbers.

            >>> from mpmath import *
            >>> mp.dps = 15
            >>> fprod([1, 2, 0.5, 7])
            mpf('7.0')

        """
        orig = ctx.prec
        try:
            v = ctx.one
            for p in factors:
                v *= p
        finally:
            ctx.prec = orig
        return +v

    def rand(ctx):
        """
        Returns an ``mpf`` with value chosen randomly from `[0, 1)`.
        The number of randomly generated bits in the mantissa is equal
        to the working precision.
        """
        return ctx.make_mpf(mpf_rand(ctx._prec))

    def fraction(ctx, p, q):
        """
        Given Python integers `(p, q)`, returns a lazy ``mpf`` representing
        the fraction `p/q`. The value is updated with the precision.

            >>> mp.dps = 15
            >>> a = fraction(1,100)
            >>> b = mpf(1)/100
            >>> print a; print b
            0.01
            0.01
            >>> mp.dps = 30
            >>> print a; print b      # a will be accurate
            0.01
            0.0100000000000000002081668171172
            >>> mp.dps = 15

        """
        return ctx.constant(lambda prec, rnd: from_rational(p, q, prec, rnd),
            '%s/%s' % (p, q))

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
            >>> mp.dps = 15
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
            >>> mp.dps = 15
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

    def mpi_from_str(ctx, s):
        """
        Parse an interval number given as a string.

        Allowed forms are
            1. 'a +- b'
            2. 'a (b%)'  % sign is optional
            3. '[a, b]'
            4. 'x[y,z]e'
        In 1, a is the midpoint of the interval and b is the half-width.
        In 2, a is the midpoint of the interval and b is the half-width.
        In 3, the interval is indicated directly.
        In 4, x are shared digits, y and z are unequal digits, e is the exponent.
        """
        e = ValueError("Improperly formed interval number '%s'" %s)
        s = s.replace(" ", "")
        if "+-" in s:
            # case 1
            n = [ctx.mpf(strip(i)) for i in s.split("+-")]
            return ctx.mpi(n[0] - n[1], n[0] + n[1])
        elif "(" in s:
            # case 2
            if s[0] == "(":  # Don't confuse with a complex number (x,y)
                return None
            if ")" not in s:
                raise e
            s = s.replace(")", "")
            percent = False
            if "%" in s:
                if s[-1] != "%":
                    raise e
                percent = True
                s = s.replace("%", "")
            a, p = [ctx.mpf(strip(i)) for i in s.split("(")]
            d = p
            if percent:
                d = a*p / 100
            return ctx.mpi(a - d, a + d)
        elif "," in s:
            if ('[' not in s) or (']' not in s):
                raise e
            if s[0] == '[':
                # case 3
                s = s.replace("[", "")
                s = s.replace("]", "")
                n = [ctx.mpf(strip(i)) for i in s.split(",")]
                return ctx.mpi(n[0], n[1])
            else:
                # case 4
                x, y = s.split('[')
                y, z = y.split(',')
                if 'e' in s:
                    z, e = z.split(']')
                else:
                    z, e = z.rstrip(']'), ''
                return ctx.mpi(x + y + e, x + z + e)
        else:
            return None

    def mpi_to_str(ctx, x, dps=None, use_spaces=True, brackets=('[', ']'),
                   mode='brackets', error_dps=4, **kwargs):
        """
        Convert a mpi interval to a string.

        **Arguments**

        *dps*
            decimal places to use for printing
        *use_spaces*
            use spaces for more readable output, defaults to true
        *brackets*
            tuple of two strings indicating the brackets to use
        *mode*
            mode of display: 'plusminus', 'percent', 'brackets' (default) or 'diff'
        *error_dps*
            limit the error to *error_dps* digits (mode 'plusminus and 'percent')

        **Examples**

            >>> from mpmath import mpi, mp
            >>> mp.dps = 30
            >>> x = mpi(1, 2)
            >>> mpi_to_str(x, mode='plusminus')
            '1.5 +- 5.0e-1'
            >>> mpi_to_str(x, mode='percent')
            '1.5 (33.33%)'
            >>> mpi_to_str(x, mode='brackets')
            '[1.0, 2.0]'
            >>> mpi_to_str(x, mode='brackets' , brackets=('<', '>'))
            '<1.0, 2.0>'
            >>> x = mpi('5.2582327113062393041', '5.2582327113062749951')
            >>> mpi_to_str(x, mode='diff')
            '5.2582327113062[4, 7]'
            >>> mpi_to_str(mpi(0), mode='percent')
            '0.0 (0%)'

        """
        if dps is None:
            dps = ctx.dps # TODO: maybe choose a smaller default value
        a = to_str(x.a._mpf_, dps, **kwargs)
        b = to_str(x.b._mpf_, dps, **kwargs)
        mid = to_str(x.mid._mpf_, dps, **kwargs)
        delta = to_str((x.delta/2)._mpf_, error_dps, **kwargs)
        sp = ""
        if use_spaces:
            sp = " "
        br1, br2 = brackets
        if mode == 'plusminus':
            s = mid + sp + "+-" + sp + delta
        elif mode == 'percent':
            a = x.mid
            if x.mid != 0:
                b = 100*x.delta/(2*x.mid)
            else:
                b = MP_ZERO
            m = str(a)
            p = ctx.nstr(b, error_dps)
            s = m + sp + "(" + p + "%)"
        elif mode == 'brackets':
            s = br1 + a.strip() + "," + sp + b + br2
        elif mode == 'diff':
            # use more digits if str(x.a) and str(x.b) are equal
            if a == b:
                a = to_str(x.a._mpf_, repr_dps(ctx.prec), **kwargs)
                b = to_str(x.b._mpf_, repr_dps(ctx.prec), **kwargs)
            # separate mantissa and exponent
            a = a.split('e')
            if len(a) == 1:
                a.append('')
            b = b.split('e')
            if len(b) == 1:
                b.append('')
            if a[1] == b[1]:
                if a[0] != b[0]:
                    for i in xrange(len(a[0]) + 1):
                        if a[0][i] != b[0][i]:
                            break
                    s = (a[0][:i] + br1 + a[0][i:] + ',' + sp + b[0][i:] + br2
                         + 'e'*min(len(a[1]), 1) + a[1])
                else: # no difference
                    s = a[0] + br1 + br2 + 'e'*min(len(a[1]), 1) + a[1]
            else:
                s = br1 + 'e'.join(a) + ',' + sp + 'e'.join(b) + br2
        else:
            raise ValueError("'%s' is unknown mode for printing mpi" % mode)
        return s



class mpnumeric(object):
    """Base class for mpf and mpc."""
    __slots__ = []
    def __new__(cls, val):
        raise NotImplementedError

get_complex = re.compile(r'^\(?(?P<re>[\+\-]?\d*\.?\d*(e[\+\-]?\d+)?)??'
                         r'(?P<im>[\+\-]?\d*\.?\d*(e[\+\-]?\d+)?j)?\)?$')
# TODO: add tests


class _mpf(mpnumeric):
    """
    An mpf instance holds a real-valued floating-point number. mpf:s
    work analogously to Python floats, but support arbitrary-precision
    arithmetic.
    """
    __slots__ = ['_mpf_']

    def __new__(cls, val=fzero, **kwargs):
        """A new mpf can be created from a Python float, an int, a
        or a decimal string representing a number in floating-point
        format."""
        prec, rounding = cls.context._prec_rounding
        if kwargs:
            prec = kwargs.get('prec', prec)
            if 'dps' in kwargs:
                prec = dps_to_prec(kwargs['dps'])
            rounding = kwargs.get('rounding', rounding)
        if type(val) is cls:
            sign, man, exp, bc = val._mpf_
            if (not man) and exp:
                return val
            v = new(cls)
            v._mpf_ = normalize(sign, man, exp, bc, prec, rounding)
            return v
        elif type(val) is tuple:
            if len(val) == 2:
                v = new(cls)
                v._mpf_ = from_man_exp(val[0], val[1], prec, rounding)
                return v
            if len(val) == 4:
                sign, man, exp, bc = val
                v = new(cls)
                v._mpf_ = normalize(sign, MP_BASE(man), exp, bc, prec, rounding)
                return v
            raise ValueError
        else:
            v = new(cls)
            v._mpf_ = mpf_pos(cls.mpf_convert_arg(val, prec, rounding), prec, rounding)
            return v

    @classmethod
    def mpf_convert_arg(cls, x, prec, rounding):
        if isinstance(x, int_types): return from_int(x)
        if isinstance(x, float): return from_float(x)
        if isinstance(x, basestring): return from_str(x, prec, rounding)
        if isinstance(x, constant): return x.func(prec, rounding)
        if hasattr(x, '_mpf_'): return x._mpf_
        if hasattr(x, '_mpmath_'):
            t = cls.context.convert(x._mpmath_(prec, rounding))
            if hasattr(t, '_mpf_'):
                return t._mpf_
        raise TypeError("cannot create mpf from " + repr(x))

    @classmethod
    def mpf_convert_rhs(cls, x):
        if isinstance(x, int_types): return from_int(x)
        if isinstance(x, float): return from_float(x)
        if isinstance(x, complex_types): return cls.context.mpc(x)
        if hasattr(x, '_mpf_'): return x._mpf_
        if hasattr(x, '_mpmath_'):
            t = cls.context.convert(x._mpmath_(*cls.context._prec_rounding))
            if hasattr(t, '_mpf_'):
                return t._mpf_
            return t
        return NotImplemented

    @classmethod
    def mpf_convert_lhs(cls, x):
        x = cls.mpf_convert_rhs(x)
        if type(x) is tuple:
            return cls.context.make_mpf(x)
        return x

    man_exp = property(lambda self: self._mpf_[1:3])
    man = property(lambda self: self._mpf_[1])
    exp = property(lambda self: self._mpf_[2])
    bc = property(lambda self: self._mpf_[3])

    real = property(lambda self: self)
    imag = property(lambda self: self.context.zero)

    conjugate = lambda self: self

    def __getstate__(self): return to_pickable(self._mpf_)
    def __setstate__(self, val): self._mpf_ = from_pickable(val)

    def __repr__(s): return "mpf('%s')" % to_str(s._mpf_, s.context.repr_digits)
    def __str__(s): return to_str(s._mpf_, s.context.str_digits)
    def __hash__(s): return mpf_hash(s._mpf_)
    def __int__(s): return int(to_int(s._mpf_))
    def __long__(s): return long(to_int(s._mpf_))
    def __float__(s): return to_float(s._mpf_)
    def __complex__(s): return complex(float(s))
    def __nonzero__(s): return s._mpf_ != fzero

    def __abs__(s):
        cls, new, (prec, rounding) = s._ctxdata
        v = new(cls)
        v._mpf_ = mpf_abs(s._mpf_, prec, rounding)
        return v

    def __pos__(s):
        cls, new, (prec, rounding) = s._ctxdata
        v = new(cls)
        v._mpf_ = mpf_pos(s._mpf_, prec, rounding)
        return v

    def __neg__(s):
        cls, new, (prec, rounding) = s._ctxdata
        v = new(cls)
        v._mpf_ = mpf_neg(s._mpf_, prec, rounding)
        return v

    def _cmp(s, t, func):
        if hasattr(t, '_mpf_'):
            t = t._mpf_
        else:
            t = s.mpf_convert_rhs(t)
            if t is NotImplemented:
                return t
        return func(s._mpf_, t)

    def __cmp__(s, t): return s._cmp(t, mpf_cmp)
    def __lt__(s, t): return s._cmp(t, mpf_lt)
    def __gt__(s, t): return s._cmp(t, mpf_gt)
    def __le__(s, t): return s._cmp(t, mpf_le)
    def __ge__(s, t): return s._cmp(t, mpf_ge)

    def __ne__(s, t):
        v = s.__eq__(t)
        if v is NotImplemented:
            return v
        return not v

    def __rsub__(s, t):
        cls, new, (prec, rounding) = s._ctxdata
        if type(t) in int_types:
            v = new(cls)
            v._mpf_ = mpf_sub(from_int(t), s._mpf_, prec, rounding)
            return v
        t = s.mpf_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t - s

    def __rdiv__(s, t):
        cls, new, (prec, rounding) = s._ctxdata
        if isinstance(t, int_types):
            v = new(cls)
            v._mpf_ = mpf_rdiv_int(t, s._mpf_, prec, rounding)
            return v
        t = s.mpf_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t / s

    def __rpow__(s, t):
        t = s.mpf_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t ** s

    def __rmod__(s, t):
        t = s.mpf_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t % s

    def sqrt(s):
        return s.context.sqrt(s)

    def ae(s, t, rel_eps=None, abs_eps=None):
        return s.context.almosteq(s, t, rel_eps, abs_eps)


mpf_binary_op = """
def %NAME%(self, other):
    mpf, new, (prec, rounding) = self._ctxdata
    sval = self._mpf_
    if hasattr(other, '_mpf_'):
        tval = other._mpf_
        %WITH_MPF%
    ttype = type(other)
    if ttype in int_types:
        %WITH_INT%
    elif ttype is float:
        tval = from_float(other)
        %WITH_MPF%
    elif hasattr(other, '_mpc_'):
        tval = other._mpc_
        mpc = type(other)
        %WITH_MPC%
    elif ttype is complex:
        tval = from_float(other.real), from_float(other.imag)
        mpc = self.context.mpc
        %WITH_MPC%
    if isinstance(other, mpnumeric):
        return NotImplemented
    try:
        other = mpmathify(other, strings=False)
    except TypeError:
        return NotImplemented
    return self.%NAME%(other)
"""

return_mpf = "; obj = new(mpf); obj._mpf_ = val; return obj"
return_mpc = "; obj = new(mpc); obj._mpc_ = val; return obj"

mpf_pow_same = """
        try:
            val = mpf_pow(sval, tval, prec, rounding) %s
        except ComplexResult:
            if mpf.context.trap_complex:
                raise
            mpc = mpf.context.mpc
            val = mpc_pow((sval, fzero), (tval, fzero), prec, rounding) %s
""" % (return_mpf, return_mpc)

def binary_op(name, with_mpf='', with_int='', with_mpc=''):
    code = mpf_binary_op
    code = code.replace("%WITH_INT%", with_int)
    code = code.replace("%WITH_MPC%", with_mpc)
    code = code.replace("%WITH_MPF%", with_mpf)
    code = code.replace("%NAME%", name)
    np = {}
    exec code in globals(), np
    return np[name]

_mpf.__eq__ = binary_op('__eq__',
    'return mpf_eq(sval, tval)',
    'return mpf_eq(sval, from_int(other))',
    'return (tval[1] == fzero) and mpf_eq(tval[0], sval)')

_mpf.__add__ = binary_op('__add__',
    'val = mpf_add(sval, tval, prec, rounding)' + return_mpf,
    'val = mpf_add(sval, from_int(other), prec, rounding)' + return_mpf,
    'val = mpc_add_mpf(tval, sval, prec, rounding)' + return_mpc)

_mpf.__sub__ = binary_op('__sub__',
    'val = mpf_sub(sval, tval, prec, rounding)' + return_mpf,
    'val = mpf_sub(sval, from_int(other), prec, rounding)' + return_mpf,
    'val = mpc_sub((sval, fzero), tval, prec, rounding)' + return_mpc)

_mpf.__mul__ = binary_op('__mul__',
    'val = mpf_mul(sval, tval, prec, rounding)' + return_mpf,
    'val = mpf_mul_int(sval, other, prec, rounding)' + return_mpf,
    'val = mpc_mul_mpf(tval, sval, prec, rounding)' + return_mpc)

_mpf.__div__ = binary_op('__div__',
    'val = mpf_div(sval, tval, prec, rounding)' + return_mpf,
    'val = mpf_div(sval, from_int(other), prec, rounding)' + return_mpf,
    'val = mpc_mpf_div(sval, tval, prec, rounding)' + return_mpc)

_mpf.__mod__ = binary_op('__mod__',
    'val = mpf_mod(sval, tval, prec, rounding)' + return_mpf,
    'val = mpf_mod(sval, from_int(other), prec, rounding)' + return_mpf,
    'raise NotImplementedError("complex modulo")')

_mpf.__pow__ = binary_op('__pow__',
    mpf_pow_same,
    'val = mpf_pow_int(sval, other, prec, rounding)' + return_mpf,
    'val = mpc_pow((sval, fzero), tval, prec, rounding)' + return_mpc)

_mpf.__radd__ = _mpf.__add__
_mpf.__rmul__ = _mpf.__mul__
_mpf.__truediv__ = _mpf.__div__
_mpf.__rtruediv__ = _mpf.__rdiv__


class _mpc(mpnumeric):
    """
    An mpc represents a complex number using a pair of mpf:s (one
    for the real part and another for the imaginary part.) The mpc
    class behaves fairly similarly to Python's complex type.
    """

    __slots__ = ['_mpc_']

    def __new__(cls, real=0, imag=0):
        s = object.__new__(cls)
        if isinstance(real, complex_types):
            real, imag = real.real, real.imag
        elif hasattr(real, '_mpc_'):
            s._mpc_ = real._mpc_
            return s
        real = cls.context.mpf(real)
        imag = cls.context.mpf(imag)
        s._mpc_ = (real._mpf_, imag._mpf_)
        return s

    real = property(lambda self: self.context.make_mpf(self._mpc_[0]))
    imag = property(lambda self: self.context.make_mpf(self._mpc_[1]))

    def __getstate__(self):
        return to_pickable(self._mpc_[0]), to_pickable(self._mpc_[1])

    def __setstate__(self, val):
        self._mpc_ = from_pickable(val[0]), from_pickable(val[1])

    def __repr__(s):
        r = repr(s.real)[4:-1]
        i = repr(s.imag)[4:-1]
        return "%s(real=%s, imag=%s)" % (type(s).__name__, r, i)

    def __str__(s):
        return "(%s)" % mpc_to_str(s._mpc_, s.context.dps)

    def __complex__(s):
        return mpc_to_complex(s._mpc_)

    def __pos__(s):
        cls, new, (prec, rounding) = s._ctxdata
        v = new(cls)
        v._mpc_ = mpc_pos(s._mpc_, prec, rounding)
        return v

    def __abs__(s):
        prec, rounding = s.context._prec_rounding
        v = new(s.context.mpf)
        v._mpf_ = mpc_abs(s._mpc_, prec, rounding)
        return v

    def __neg__(s):
        cls, new, (prec, rounding) = s._ctxdata
        v = new(cls)
        v._mpc_ = mpc_neg(s._mpc_, prec, rounding)
        return v

    def conjugate(s):
        cls, new, (prec, rounding) = s._ctxdata
        v = new(cls)
        v._mpc_ = mpc_conjugate(s._mpc_, prec, rounding)
        return v

    def __nonzero__(s):
        return mpc_is_nonzero(s._mpc_)

    def __hash__(s):
        return mpc_hash(s._mpc_)

    @classmethod
    def mpc_convert_lhs(cls, x):
        try:
            y = cls.context.convert(x)
            return y
        except TypeError:
            return NotImplemented

    def __eq__(s, t):
        if not hasattr(t, '_mpc_'):
            if isinstance(t, str):
                return False
            t = s.mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
        return s.real == t.real and s.imag == t.imag

    def __ne__(s, t):
        b = s.__eq__(t)
        if b is NotImplemented:
            return b
        return not b

    def _compare(*args):
        raise TypeError("no ordering relation is defined for complex numbers")

    __gt__ = _compare
    __le__ = _compare
    __gt__ = _compare
    __ge__ = _compare

    def __add__(s, t):
        cls, new, (prec, rounding) = s._ctxdata
        if not hasattr(t, '_mpc_'):
            t = s.mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
            if hasattr(t, '_mpf_'):
                v = new(cls)
                v._mpc_ = mpc_add_mpf(s._mpc_, t._mpf_, prec, rounding)
                return v
        v = new(cls)
        v._mpc_ = mpc_add(s._mpc_, t._mpc_, prec, rounding)
        return v

    def __sub__(s, t):
        cls, new, (prec, rounding) = s._ctxdata
        if not hasattr(t, '_mpc_'):
            t = s.mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
            if hasattr(t, '_mpf_'):
                v = new(cls)
                v._mpc_ = mpc_sub_mpf(s._mpc_, t._mpf_, prec, rounding)
                return v
        v = new(cls)
        v._mpc_ = mpc_sub(s._mpc_, t._mpc_, prec, rounding)
        return v

    def __mul__(s, t):
        cls, new, (prec, rounding) = s._ctxdata
        if not hasattr(t, '_mpc_'):
            if isinstance(t, int_types):
                v = new(cls)
                v._mpc_ = mpc_mul_int(s._mpc_, t, prec, rounding)
                return v
            t = s.mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
            if hasattr(t, '_mpf_'):
                v = new(cls)
                v._mpc_ = mpc_mul_mpf(s._mpc_, t._mpf_, prec, rounding)
            t = mpc(t)
        v = new(cls)
        v._mpc_ = mpc_mul(s._mpc_, t._mpc_, prec, rounding)
        return v

    def __div__(s, t):
        cls, new, (prec, rounding) = s._ctxdata
        if not hasattr(t, '_mpc_'):
            t = s.mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
            if hasattr(t, '_mpf_'):
                v = new(cls)
                v._mpc_ = mpc_div_mpf(s._mpc_, t._mpf_, prec, rounding)
                return v
        v = new(cls)
        v._mpc_ = mpc_div(s._mpc_, t._mpc_, prec, rounding)
        return v

    def __pow__(s, t):
        cls, new, (prec, rounding) = s._ctxdata
        if isinstance(t, int_types):
            v = new(cls)
            v._mpc_ = mpc_pow_int(s._mpc_, t, prec, rounding)
            return v
        t = s.mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        v = new(cls)
        if hasattr(t, '_mpf_'):
            v._mpc_ = mpc_pow_mpf(s._mpc_, t._mpf_, prec, rounding)
        else:
            v._mpc_ = mpc_pow(s._mpc_, t._mpc_, prec, rounding)
        return v

    __radd__ = __add__

    def __rsub__(s, t):
        t = s.mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t - s

    def __rmul__(s, t):
        cls, new, (prec, rounding) = s._ctxdata
        if isinstance(t, int_types):
            v = new(cls)
            v._mpc_ = mpc_mul_int(s._mpc_, t, prec, rounding)
            return v
        t = s.mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t * s

    def __rdiv__(s, t):
        t = s.mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t / s

    def __rpow__(s, t):
        t = s.mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t ** s

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def ae(s, t, rel_eps=None, abs_eps=None):
        return s.context.almosteq(s, t, rel_eps, abs_eps)


class _mpi(mpnumeric):
    """
    Interval arithmetic class. Precision is controlled by mp.prec.
    """

    def __new__(cls, a, b=None):
        if isinstance(a, mpi):
            return a
        if b is None:
            b = a
        a = mpf(a, rounding=round_floor)
        b = mpf(b, rounding=round_ceiling)
        if isnan(a) or isnan(b):
            a, b = -inf, inf
        assert a <= b, "endpoints must be properly ordered"
        return make_mpi((a._mpf_, b._mpf_))

    @property
    def a(self):
        return make_mpf(self._mpi_[0])

    @property
    def b(self):
        return make_mpf(self._mpi_[1])

    @property
    def mid(self):
        return make_mpf(mpi_mid(self._mpi_, mp.prec))

    @property
    def delta(self):
        return make_mpf(mpi_delta(self._mpi_, mp.prec))

    def _compare(*args):
        raise TypeError("no ordering relation is defined for intervals")

    __gt__ = _compare
    __le__ = _compare
    __gt__ = _compare
    __ge__ = _compare

    def __contains__(self, t):
        t = mpi(t)
        return (self.a <= t.a) and (t.b <= self.b)

    def __str__(self):
        return mpi_str(self._mpi_, mp.prec)

    def __repr__(self):
        return "mpi(%r, %r)" % (self.a, self.b)

    def __eq__(self, other):
        if not isinstance(other, mpi):
            try:
                other = mpi(other)
            except:
                return NotImplemented
        return (self.a == other.a) and (self.b == other.b)

    def __ne__(self, other):
        return not (self == other)

    def __abs__(self):
        return make_mpi(mpi_abs(self._mpi_, mp.prec))

    def __pos__(self):
        return make_mpi(mpi_pos(self._mpi_, mp.prec))

    def __neg__(self):
        return make_mpi(mpi_neg(self._mpi_, mp.prec))

    def __add__(self, other):
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_add(self._mpi_, other._mpi_, mp.prec))

    def __sub__(self, other):
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_sub(self._mpi_, other._mpi_, mp.prec))

    def __mul__(self, other):
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_mul(self._mpi_, other._mpi_, mp.prec))

    def __div__(self, other):
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_div(self._mpi_, other._mpi_, mp.prec))

    def __pow__(self, other):
        if isinstance(other, (int, long)):
            return make_mpi(mpi_pow_int(self._mpi_, int(other), mp.prec))
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_pow(self._mpi_, other._mpi_, mp.prec))

    def __rsub__(s, t):
        return mpi(t) - s

    def __rdiv__(s, t):
        return mpi(t) / s

    def __rpow__(s, t):
        return mpi(t) ** s

    __radd__ = __add__
    __rmul__ = __mul__
    __truediv__ = __div__
    __rtruediv__ = __rdiv__
    __floordiv__ = __div__
    __rfloordiv__ = __rdiv__



class _constant(_mpf):
    """Represents a mathematical constant with dynamic precision.
    When printed or used in an arithmetic operation, a constant
    is converted to a regular mpf at the working precision. A
    regular mpf can also be obtained using the operation +x."""

    def __new__(cls, func, name):
        a = object.__new__(cls)
        a.name = name
        a.func = func
        return a

    def __call__(self, prec=None, dps=None, rounding=None):
        prec2, rounding2 = self.context._prec_rounding
        if not prec: prec = prec2
        if not rounding: rounding = rounding2
        if dps: prec = dps_to_prec(dps)
        return self.context.make_mpf(self.func(prec, rounding))

    @property
    def _mpf_(self):
        prec, rounding = self.context._prec_rounding
        return self.func(prec, rounding)

    def __repr__(self):
        return "<%s: %s~>" % (self.name, self.context.nstr(self))


def absmin(x):
    """
    Returns ``abs(x).a`` for an interval, or ``abs(x)`` for anything else.
    """
    if hasattr(x, '_mpi_'):
        return abs(x).a
    return abs(x)

def absmax(x):
    """
    Returns ``abs(x).b`` for an interval, or ``abs(x)`` for anything else.
    """
    if hasattr(x, '_mpi_'):
        return abs(x).b
    return abs(x)

def AS_POINTS(x):
    if hasattr(x, '_mpi_'):
        return [x.a, x.b]
    return x


def def_mp_builtin(name, mpf_f, mpc_f=None, mpi_f=None, doc="<no doc>"):
    """
    Given a low-level mpf_ function, and optionally similar functions
    for mpc_ and mpi_, defines the function as a MultiPrecisionArithmetic
    method.

    It is assumed that the return type is the same as that of
    the input; the exception is that propagation from mpf to mpc is possible
    by raising ComplexResult.

    """
    def f(ctx, x, **kwargs):
        if type(x) not in ctx.types:
            x = ctx.convert(x)
        prec, rounding = ctx._prec_rounding
        if kwargs:
            prec = kwargs.get('prec', prec)
            if 'dps' in kwargs:
                prec = dps_to_prec(kwargs['dps'])
            rounding = kwargs.get('rounding', rounding)
        if hasattr(x, '_mpf_'):
            try:
                return ctx.make_mpf(mpf_f(x._mpf_, prec, rounding))
            except ComplexResult:
                # Handle propagation to complex
                if ctx.trap_complex:
                    raise
                return ctx.make_mpc(mpc_f((x._mpf_, fzero), prec, rounding))
        elif hasattr(x, '_mpc_'):
            return ctx.make_mpc(mpc_f(x._mpc_, prec, rounding))
        elif hasattr(x, '_mpi_'):
            if mpi_f:
                return ctx.make_mpi(mpi_f(x._mpi_, prec))
        raise NotImplementedError("%s of a %s" % (name, type(x)))
    f.__name__ = name
    import function_docs
    f.__doc__ = function_docs.__dict__.get(name, "Computes the %s of x" % doc)
    setattr(MultiPrecisionArithmetic, name, f)
    return getattr(MultiPrecisionArithmetic.default_instance, name)

def defun_wrapped(f):
    """
    Defines function as a MultiPrecisionArithmetic method, and adds
    wrapper code for setting precision (including allocating some
    extra working precision before entry) as well as converting
    all input args to recognized number types.
    """
    def g(ctx, *args, **kwargs):
        orig_prec = ctx.prec
        try:
            ctx.prec = orig_prec + 10
            args = [ctx.convert(z) for z in args]
            v = f(ctx, *args, **kwargs)
        finally:
            ctx.prec = orig_prec
        return +v
    g.__name__ = f.__name__
    import function_docs
    g.__doc__ = function_docs.__dict__.get(f.__name__, f.__doc__)
    setattr(MultiPrecisionArithmetic, f.__name__, g)
    return getattr(MultiPrecisionArithmetic.default_instance, f.__name__)

def defun(f):
    """
    Defines function as a MultiPrecisionArithmetic method.
    """
    import function_docs
    f.__doc__ = function_docs.__dict__.get(f.__name__, f.__doc__)
    setattr(MultiPrecisionArithmetic, f.__name__, f)
    return getattr(MultiPrecisionArithmetic.default_instance, f.__name__)

def defun_static(f):
    """
    Defines function as a MultiPrecisionArithmetic static method.
    """
    import function_docs
    f.__doc__ = function_docs.__dict__.get(f.__name__, f.__doc__)
    setattr(MultiPrecisionArithmetic, f.__name__, staticmethod(f))
    return f

complex_types = (complex, _mpc)


# Standard instance
mp = MultiPrecisionArithmetic()
mp.default()

MultiPrecisionArithmetic.default_instance = mp

make_mpi = mp.make_mpi
make_mpf = mp.make_mpf
make_mpc = mp.make_mpc

prec_rounding = mp._prec_rounding

workprec = mp.workprec
workdps = mp.workdps
extraprec = mp.extraprec
extradps = mp.extradps

eps = mp.eps
one = mp.one
zero = mp.zero
inf = mp.inf
ninf = mp.ninf
nan = mp.nan
j = mp.j
isnan = mp.isnan
isinf = mp.isinf
isint = mp.isint
almosteq = mp.almosteq
chop = mp.chop
mpf = mp.mpf
mpc = mp.mpc
mpi = mp.mpi
constant = mp.constant
mpmathify = mp.convert
nstr = mp.nstr
nprint = mp.nprint
rand = mp.rand
fsum = mp.fsum
fdot = mp.fdot
fprod = mp.fprod
fraction = mp.fraction
arange = mp.arange
linspace = mp.linspace

mpi_to_str = mp.mpi_to_str
mpi_from_str = mp.mpi_from_str

# XXX
mp.absmin = absmin
mp.absmax = absmax

if __name__ == '__main__':
    import doctest
    doctest.testmod()
