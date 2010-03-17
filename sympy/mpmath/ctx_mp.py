"""
This module defines the mpf, mpc classes, and standard functions for
operating with them.
"""
__docformat__ = 'plaintext'

import re

from string import strip

from ctx_base import StandardBaseContext

import libmp

from libmp import (MPZ, MPZ_ZERO, MPZ_ONE, int_types, repr_dps,
    round_floor, round_ceiling, dps_to_prec, round_nearest, prec_to_dps,
    ComplexResult, to_pickable, from_pickable, normalize,
    from_int, from_float, from_str, to_int, to_float, to_str,
    from_rational, from_man_exp,
    fone, fzero, finf, fninf, fnan,
    mpf_abs, mpf_pos, mpf_neg, mpf_add, mpf_sub, mpf_mul, mpf_mul_int,
    mpf_div, mpf_rdiv_int, mpf_pow_int, mpf_mod,
    mpf_eq, mpf_cmp, mpf_lt, mpf_gt, mpf_le, mpf_ge,
    mpf_hash, mpf_rand,
    mpf_sum,
    bitcount, to_fixed,
    mpc_to_str,
    mpc_to_complex, mpc_hash, mpc_pos, mpc_is_nonzero, mpc_neg, mpc_conjugate,
    mpc_abs, mpc_add, mpc_add_mpf, mpc_sub, mpc_sub_mpf, mpc_mul, mpc_mul_mpf,
    mpc_mul_int, mpc_div, mpc_div_mpf, mpc_pow, mpc_pow_mpf, mpc_pow_int,
    mpc_mpf_div,
    mpf_pow,
    mpi_mid, mpi_delta, mpi_str,
    mpi_abs, mpi_pos, mpi_neg, mpi_add, mpi_sub,
    mpi_mul, mpi_div, mpi_pow_int, mpi_pow,
    mpf_pi, mpf_degree, mpf_e, mpf_phi, mpf_ln2, mpf_ln10,
    mpf_euler, mpf_catalan, mpf_apery, mpf_khinchin,
    mpf_glaisher, mpf_twinprime, mpf_mertens,
    int_types)

import function_docs
import rational

new = object.__new__

get_complex = re.compile(r'^\(?(?P<re>[\+\-]?\d*\.?\d*(e[\+\-]?\d+)?)??'
                         r'(?P<im>[\+\-]?\d*\.?\d*(e[\+\-]?\d+)?j)?\)?$')


try:
    from sage.libs.mpmath.ext_main import Context as BaseMPContext
    # pickle hack
    import sage.libs.mpmath.ext_main as _mpf_module
except ImportError:
    from ctx_mp_python import PythonMPContext as BaseMPContext
    import ctx_mp_python as _mpf_module

from ctx_mp_python import _mpf, _mpc, mpnumeric


class _mpi(mpnumeric):
    """
    Interval arithmetic class. Precision is controlled by mp.prec.
    """

    def __new__(cls, a, b=None):
        ctx = cls.context
        if isinstance(a, ctx.mpi):
            return a
        if b is None:
            b = a
        a = ctx.mpf(a, rounding=round_floor)
        b = ctx.mpf(b, rounding=round_ceiling)
        if ctx.isnan(a) or ctx.isnan(b):
            a, b = ctx.ninf, ctx.inf
        assert a <= b, "endpoints must be properly ordered"
        return ctx.make_mpi((a._mpf_, b._mpf_))

    @property
    def a(self):
        return self.context.make_mpf(self._mpi_[0])

    @property
    def b(self):
        return self.context.make_mpf(self._mpi_[1])

    @property
    def mid(self):
        ctx = self.context
        return ctx.make_mpf(mpi_mid(self._mpi_, ctx.prec))

    @property
    def delta(self):
        ctx = self.context
        return ctx.make_mpf(mpi_delta(self._mpi_, ctx.prec))

    def _compare(*args):
        raise TypeError("no ordering relation is defined for intervals")

    __gt__ = _compare
    __le__ = _compare
    __gt__ = _compare
    __ge__ = _compare

    def __contains__(self, t):
        t = self.context.mpi(t)
        return (self.a <= t.a) and (t.b <= self.b)

    def __str__(self):
        return mpi_str(self._mpi_, self.context.prec)

    def __repr__(self):
        if self.context.pretty:
            return str(self)
        return "mpi(%r, %r)" % (self.a, self.b)

    def __eq__(self, other):
        if not hasattr(other, "_mpi_"):
            try:
                other = self.context.mpi(other)
            except:
                return NotImplemented
        return (self.a == other.a) and (self.b == other.b)

    def __ne__(self, other):
        return not (self == other)

    def __abs__(self):
        return self.context.make_mpi(mpi_abs(self._mpi_, self.context.prec))

    def __pos__(self):
        return self.context.make_mpi(mpi_pos(self._mpi_, self.context.prec))

    def __neg__(self):
        return self.context.make_mpi(mpi_neg(self._mpi_, self.context.prec))

    def __add__(self, other):
        if not hasattr(other, "_mpi_"):
            other = self.context.mpi(other)
        return self.context.make_mpi(mpi_add(self._mpi_, other._mpi_,
            self.context.prec))

    def __sub__(self, other):
        if not hasattr(other, "_mpi_"):
            other = self.context.mpi(other)
        return self.context.make_mpi(mpi_sub(self._mpi_, other._mpi_,
            self.context.prec))

    def __mul__(self, other):
        if not hasattr(other, "_mpi_"):
            other = self.context.mpi(other)
        return self.context.make_mpi(mpi_mul(self._mpi_, other._mpi_,
            self.context.prec))

    def __div__(self, other):
        if not hasattr(other, "_mpi_"):
            other = self.context.mpi(other)
        return self.context.make_mpi(mpi_div(self._mpi_, other._mpi_,
            self.context.prec))

    def __pow__(self, other):
        if isinstance(other, (int, long)):
            return self.context.make_mpi(mpi_pow_int(self._mpi_, int(other),
                self.context.prec))
        if not hasattr(other, "_mpi_"):
            other = self.context.mpi(other)
        return self.context.make_mpi(mpi_pow(self._mpi_, other._mpi_,
            self.context.prec))

    def __rsub__(s, t):
        return s.context.mpi(t) - s

    def __rdiv__(s, t):
        return s.context.mpi(t) / s

    def __rpow__(s, t):
        return s.context.mpi(t) ** s

    __radd__ = __add__
    __rmul__ = __mul__
    __truediv__ = __div__
    __rtruediv__ = __rdiv__
    __floordiv__ = __div__
    __rfloordiv__ = __rdiv__

class MPContext(BaseMPContext, StandardBaseContext):
    """
    Context for multiprecision arithmetic with a global precision.
    """

    def __init__(ctx):
        BaseMPContext.__init__(ctx)

        ctx.trap_complex = False
        ctx.pretty = False
        ctx.mpi = type('mpi', (_mpi,), {})
        ctx.types = [ctx.mpf, ctx.mpc, ctx.mpi, ctx.constant]
        # For fast access
        ctx.mpi._ctxdata = [ctx.mpi, new, ctx._prec_rounding]
        ctx.mpi.context = ctx

        ctx._mpq = rational.mpq

        ctx.default()
        StandardBaseContext.__init__(ctx)

        ctx.mpq = rational.mpq
        ctx.init_builtins()

        ctx.hyp_summators = {}

        ctx._init_aliases()

        # XXX: automate
        ctx.bernoulli.im_func.func_doc = function_docs.bernoulli
        ctx.primepi.im_func.func_doc = function_docs.primepi
        ctx.psi.im_func.func_doc = function_docs.psi
        ctx.atan2.im_func.func_doc = function_docs.atan2
        ctx.digamma.func_doc = function_docs.digamma
        ctx.cospi.func_doc = function_docs.cospi
        ctx.sinpi.func_doc = function_docs.sinpi

    def init_builtins(ctx):

        mpf = ctx.mpf
        mpc = ctx.mpc

        # Exact constants
        ctx.one = ctx.make_mpf(fone)
        ctx.zero = ctx.make_mpf(fzero)
        ctx.j = ctx.make_mpc((fzero,fone))
        ctx.inf = ctx.make_mpf(finf)
        ctx.ninf = ctx.make_mpf(fninf)
        ctx.nan = ctx.make_mpf(fnan)

        eps = ctx.constant(lambda prec, rnd: (0, MPZ_ONE, 1-prec, 1),
            "epsilon of working precision", "eps")
        ctx.eps = eps

        # Approximate constants
        ctx.pi = ctx.constant(mpf_pi, "pi", "pi")
        ctx.ln2 = ctx.constant(mpf_ln2, "ln(2)", "ln2")
        ctx.ln10 = ctx.constant(mpf_ln10, "ln(10)", "ln10")
        ctx.phi = ctx.constant(mpf_phi, "Golden ratio phi", "phi")
        ctx.e = ctx.constant(mpf_e, "e = exp(1)", "e")
        ctx.euler = ctx.constant(mpf_euler, "Euler's constant", "euler")
        ctx.catalan = ctx.constant(mpf_catalan, "Catalan's constant", "catalan")
        ctx.khinchin = ctx.constant(mpf_khinchin, "Khinchin's constant", "khinchin")
        ctx.glaisher = ctx.constant(mpf_glaisher, "Glaisher's constant", "glaisher")
        ctx.apery = ctx.constant(mpf_apery, "Apery's constant", "apery")
        ctx.degree = ctx.constant(mpf_degree, "1 deg = pi / 180", "degree")
        ctx.twinprime = ctx.constant(mpf_twinprime, "Twin prime constant", "twinprime")
        ctx.mertens = ctx.constant(mpf_mertens, "Mertens' constant", "mertens")

        # Standard functions
        ctx.sqrt = ctx._wrap_libmp_function(libmp.mpf_sqrt, libmp.mpc_sqrt, libmp.mpi_sqrt)
        ctx.cbrt = ctx._wrap_libmp_function(libmp.mpf_cbrt, libmp.mpc_cbrt)
        ctx.ln = ctx._wrap_libmp_function(libmp.mpf_log, libmp.mpc_log, libmp.mpi_log)
        ctx.atan = ctx._wrap_libmp_function(libmp.mpf_atan, libmp.mpc_atan)
        ctx.exp = ctx._wrap_libmp_function(libmp.mpf_exp, libmp.mpc_exp, libmp.mpi_exp)
        ctx.expj = ctx._wrap_libmp_function(libmp.mpf_expj, libmp.mpc_expj)
        ctx.expjpi = ctx._wrap_libmp_function(libmp.mpf_expjpi, libmp.mpc_expjpi)
        ctx.sin = ctx._wrap_libmp_function(libmp.mpf_sin, libmp.mpc_sin, libmp.mpi_sin)
        ctx.cos = ctx._wrap_libmp_function(libmp.mpf_cos, libmp.mpc_cos, libmp.mpi_cos)
        ctx.tan = ctx._wrap_libmp_function(libmp.mpf_tan, libmp.mpc_tan, libmp.mpi_tan)
        ctx.sinh = ctx._wrap_libmp_function(libmp.mpf_sinh, libmp.mpc_sinh)
        ctx.cosh = ctx._wrap_libmp_function(libmp.mpf_cosh, libmp.mpc_cosh)
        ctx.tanh = ctx._wrap_libmp_function(libmp.mpf_tanh, libmp.mpc_tanh)
        ctx.asin = ctx._wrap_libmp_function(libmp.mpf_asin, libmp.mpc_asin)
        ctx.acos = ctx._wrap_libmp_function(libmp.mpf_acos, libmp.mpc_acos)
        ctx.atan = ctx._wrap_libmp_function(libmp.mpf_atan, libmp.mpc_atan)
        ctx.asinh = ctx._wrap_libmp_function(libmp.mpf_asinh, libmp.mpc_asinh)
        ctx.acosh = ctx._wrap_libmp_function(libmp.mpf_acosh, libmp.mpc_acosh)
        ctx.atanh = ctx._wrap_libmp_function(libmp.mpf_atanh, libmp.mpc_atanh)
        ctx.sinpi = ctx._wrap_libmp_function(libmp.mpf_sin_pi, libmp.mpc_sin_pi)
        ctx.cospi = ctx._wrap_libmp_function(libmp.mpf_cos_pi, libmp.mpc_cos_pi)
        ctx.floor = ctx._wrap_libmp_function(libmp.mpf_floor, libmp.mpc_floor)
        ctx.ceil = ctx._wrap_libmp_function(libmp.mpf_ceil, libmp.mpc_ceil)
        ctx.fib = ctx.fibonacci = ctx._wrap_libmp_function(libmp.mpf_fibonacci, libmp.mpc_fibonacci)
        ctx.gamma = ctx._wrap_libmp_function(libmp.mpf_gamma, libmp.mpc_gamma)
        ctx.digamma = ctx._wrap_libmp_function(libmp.mpf_psi0, libmp.mpc_psi0)
        ctx.fac = ctx.factorial = ctx._wrap_libmp_function(libmp.mpf_factorial, libmp.mpc_factorial)
        ctx.harmonic = ctx._wrap_libmp_function(libmp.mpf_harmonic, libmp.mpc_harmonic)
        ctx.ei = ctx._wrap_libmp_function(libmp.mpf_ei, libmp.mpc_ei)
        ctx.e1 = ctx._wrap_libmp_function(libmp.mpf_e1, libmp.mpc_e1)
        ctx._ci = ctx._wrap_libmp_function(libmp.mpf_ci, libmp.mpc_ci)
        ctx._si = ctx._wrap_libmp_function(libmp.mpf_si, libmp.mpc_si)
        ctx.ellipk = ctx._wrap_libmp_function(libmp.mpf_ellipk, libmp.mpc_ellipk)
        ctx.ellipe = ctx._wrap_libmp_function(libmp.mpf_ellipe, libmp.mpc_ellipe)
        ctx.agm1 = ctx._wrap_libmp_function(libmp.mpf_agm1, libmp.mpc_agm1)
        ctx._erf = ctx._wrap_libmp_function(libmp.mpf_erf, None)
        ctx._erfc = ctx._wrap_libmp_function(libmp.mpf_erfc, None)
        ctx._zeta = ctx._wrap_libmp_function(libmp.mpf_zeta, libmp.mpc_zeta)
        ctx._altzeta = ctx._wrap_libmp_function(libmp.mpf_altzeta, libmp.mpc_altzeta)

    def to_fixed(ctx, x, prec):
        return x.to_fixed(prec)

    def hypot(ctx, x, y):
        r"""
        Computes the Euclidean norm of the vector `(x, y)`, equal
        to `\sqrt{x^2 + y^2}`. Both `x` and `y` must be real."""
        x = ctx.convert(x)
        y = ctx.convert(y)
        return ctx.make_mpf(libmp.mpf_hypot(x._mpf_, y._mpf_, *ctx._prec_rounding))

    def _gamma_upper_int(ctx, n, z):
        n = int(n)
        if n == 0:
            return ctx.e1(z)
        if not hasattr(z, '_mpf_'):
            raise NotImplementedError
        prec, rounding = ctx._prec_rounding
        real, imag = libmp.mpf_expint(n, z._mpf_, prec, rounding, gamma=True)
        if imag is None:
            return ctx.make_mpf(real)
        else:
            return ctx.make_mpc((real, imag))

    def _expint_int(ctx, n, z):
        n = int(n)
        if n == 1:
            return ctx.e1(z)
        if not hasattr(z, '_mpf_'):
            raise NotImplementedError
        prec, rounding = ctx._prec_rounding
        real, imag = libmp.mpf_expint(n, z._mpf_, prec, rounding)
        if imag is None:
            return ctx.make_mpf(real)
        else:
            return ctx.make_mpc((real, imag))

    def _nthroot(ctx, x, n):
        if hasattr(x, '_mpf_'):
            try:
                return ctx.make_mpf(libmp.mpf_nthroot(x._mpf_, n, *ctx._prec_rounding))
            except ComplexResult:
                if ctx.trap_complex:
                    raise
                x = (x._mpf_, libmp.fzero)
        else:
            x = x._mpc_
        return ctx.make_mpc(libmp.mpc_nthroot(x, n, *ctx._prec_rounding))

    def _besselj(ctx, n, z):
        prec, rounding = ctx._prec_rounding
        if hasattr(z, '_mpf_'):
            return ctx.make_mpf(libmp.mpf_besseljn(n, z._mpf_, prec, rounding))
        elif hasattr(z, '_mpc_'):
            return ctx.make_mpc(libmp.mpc_besseljn(n, z._mpc_, prec, rounding))

    def _agm(ctx, a, b=1):
        prec, rounding = ctx._prec_rounding
        if hasattr(a, '_mpf_') and hasattr(b, '_mpf_'):
            try:
                v = libmp.mpf_agm(a._mpf_, b._mpf_, prec, rounding)
                return ctx.make_mpf(v)
            except ComplexResult:
                pass
        if hasattr(a, '_mpf_'): a = (a._mpf_, libmp.fzero)
        else: a = a._mpc_
        if hasattr(b, '_mpf_'): b = (b._mpf_, libmp.fzero)
        else: b = b._mpc_
        return ctx.make_mpc(libmp.mpc_agm(a, b, prec, rounding))

    def bernoulli(ctx, n):
        return ctx.make_mpf(libmp.mpf_bernoulli(int(n), *ctx._prec_rounding))

    def _zeta_int(ctx, n):
        return ctx.make_mpf(libmp.mpf_zeta_int(int(n), *ctx._prec_rounding))

    def atan2(ctx, y, x):
        x = ctx.convert(x)
        y = ctx.convert(y)
        return ctx.make_mpf(libmp.mpf_atan2(y._mpf_, x._mpf_, *ctx._prec_rounding))

    def psi(ctx, m, z):
        z = ctx.convert(z)
        m = int(m)
        if ctx._is_real_type(z):
            return ctx.make_mpf(libmp.mpf_psi(m, z._mpf_, *ctx._prec_rounding))
        else:
            return ctx.make_mpc(libmp.mpc_psi(m, z._mpc_, *ctx._prec_rounding))

    def clone(ctx):
        """
        Create a copy of the context, with the same working precision.
        """
        a = ctx.__class__()
        a.prec = ctx.prec
        return a

    # Several helper methods
    # TODO: add more of these, make consistent, write docstrings, ...

    def _is_real_type(ctx, x):
        if hasattr(x, '_mpc_') or type(x) is complex:
            return False
        return True

    def _is_complex_type(ctx, x):
        if hasattr(x, '_mpc_') or type(x) is complex:
            return True
        return False

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
        if type(x) in int_types:
            return x <= 0
        if isinstance(x, ctx.mpq):
            # XXX: WRONG
            p, q = x
            if not p:
                return True
            return (not (q % p)) and p <= 0
        return ctx.isnpint(ctx.convert(x))

    def __str__(ctx):
        lines = ["Mpmath settings:",
            ("  mp.prec = %s" % ctx.prec).ljust(30) + "[default: 53]",
            ("  mp.dps = %s" % ctx.dps).ljust(30) + "[default: 15]",
            ("  mp.trap_complex = %s" % ctx.trap_complex).ljust(30) + "[default: False]",
        ]
        return "\n".join(lines)

    @property
    def _repr_digits(ctx):
        return repr_dps(ctx._prec)

    @property
    def _str_digits(ctx):
        return ctx._dps

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
            return "[%s]" % (", ".join(ctx.nstr(c, n, **kwargs) for c in x))
        if isinstance(x, tuple):
            return "(%s)" % (", ".join(ctx.nstr(c, n, **kwargs) for c in x))
        if hasattr(x, '_mpf_'):
            return to_str(x._mpf_, n, **kwargs)
        if hasattr(x, '_mpc_'):
            return "(" + mpc_to_str(x._mpc_, n, **kwargs)  + ")"
        if isinstance(x, basestring):
            return repr(x)
        if isinstance(x, ctx.matrix):
            return x.__nstr__(n, **kwargs)
        if hasattr(x, '_mpi_'):
            return ctx.mpi_to_str(x, n, **kwargs)
        return str(x)

    def nprint(ctx, x, n=6, **kwargs):
        """
        Equivalent to ``print nstr(x, n)``.
        """
        print ctx.nstr(x, n, **kwargs)

    def _convert_fallback(ctx, x, strings):
        if strings and isinstance(x, basestring):
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
        if type(x) in ctx.types:  # XXX fix for mpi for Cython context
            return x
        raise TypeError("cannot create mpf from " + repr(x))

    def mpmathify(ctx, *args, **kwargs):
        return ctx.convert(*args, **kwargs)

    def _parse_prec(ctx, kwargs):
        if kwargs:
            if kwargs.get('exact'):
                return 0, 'f'
            prec, rounding = ctx._prec_rounding
            if 'rounding' in kwargs:
                rounding = kwargs['rounding']
            if 'prec' in kwargs:
                prec = kwargs['prec']
                if prec == ctx.inf:
                    return 0, 'f'
                else:
                    prec = int(prec)
            elif 'dps' in kwargs:
                dps = kwargs['dps']
                if dps == ctx.inf:
                    return 0, 'f'
                prec = dps_to_prec(dps)
            return prec, rounding
        return ctx._prec_rounding

    _exact_overflow_msg = "the exact result does not fit in memory"

    _hypsum_msg = """hypsum() failed to converge to the requested %i bits of accuracy
using a working precision of %i bits. Try with a higher maxprec,
maxterms, or set zeroprec."""

    def hypsum(ctx, p, q, flags, coeffs, z, accurate_small=True, **kwargs):
        if hasattr(z, "_mpf_"):
            key = p, q, flags, 'R'
            v = z._mpf_
        elif hasattr(z, "_mpc_"):
            key = p, q, flags, 'C'
            v = z._mpc_
        if key not in ctx.hyp_summators:
            ctx.hyp_summators[key] = libmp.make_hyp_summator(key)[1]
        summator = ctx.hyp_summators[key]
        prec = ctx.prec
        maxprec = kwargs.get('maxprec', ctx._default_hyper_maxprec(prec))
        extraprec = 50
        epsshift = 25
        # Jumps in magnitude occur when parameters are close to negative
        # integers. We must ensure that these terms are included in
        # the sum and added accurately
        magnitude_check = {}
        max_total_jump = 0
        for i, c in enumerate(coeffs):
            if flags[i] == 'Z':
                if i >= p and c <= 0:
                    ok = False
                    for ii, cc in enumerate(coeffs[:p]):
                        # Note: c <= cc or c < cc, depending on convention
                        if flags[ii] == 'Z' and cc <= 0 and c <= cc:
                            ok = True
                    if not ok:
                        raise ZeroDivisionError("pole in hypergeometric series")
                continue
            n, d = ctx.nint_distance(c)
            n = -int(n)
            d = -d
            if i >= p and n >= 0 and d > 4:
                if n in magnitude_check:
                    magnitude_check[n] += d
                else:
                    magnitude_check[n] = d
                extraprec = max(extraprec, d - prec + 60)
            max_total_jump += abs(d)
        while 1:
            if extraprec > maxprec:
                raise ValueError(ctx._hypsum_msg % (prec, prec+extraprec))
            wp = prec + extraprec
            if magnitude_check:
                mag_dict = dict((n,None) for n in magnitude_check)
            else:
                mag_dict = {}
            zv, have_complex, magnitude = summator(coeffs, v, prec, wp, \
                epsshift, mag_dict, **kwargs)
            cancel = -magnitude
            jumps_resolved = True
            if extraprec < max_total_jump:
                for n in mag_dict.values():
                    if (n is None) or (n < prec):
                        jumps_resolved = False
                        break
            accurate = (cancel < extraprec-25-5 or not accurate_small)
            if jumps_resolved:
                if accurate:
                    break
                # zero?
                zeroprec = kwargs.get('zeroprec')
                if zeroprec is not None:
                    if cancel > zeroprec:
                        if have_complex:
                            return ctx.mpc(0)
                        else:
                            return ctx.zero

            # Some near-singularities were not included, so increase
            # precision and repeat until they are
            extraprec *= 2
            # Possible workaround for bad roundoff in fixed-point arithmetic
            epsshift += 5
            extraprec += 5

        if have_complex:
            z = ctx.make_mpc(zv)
        else:
            z = ctx.make_mpf(zv)
        return z

    def ldexp(ctx, x, n):
        r"""
        Computes `x 2^n` efficiently. No rounding is performed.
        The argument `x` must be a real floating-point number (or
        possible to convert into one) and `n` must be a Python ``int``.

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> ldexp(1, 10)
            mpf('1024.0')
            >>> ldexp(1, -3)
            mpf('0.125')

        """
        x = ctx.convert(x)
        return ctx.make_mpf(libmp.mpf_shift(x._mpf_, n))

    def frexp(ctx, x):
        r"""
        Given a real number `x`, returns `(y, n)` with `y \in [0.5, 1)`,
        `n` a Python integer, and such that `x = y 2^n`. No rounding is
        performed.

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> frexp(7.5)
            (mpf('0.9375'), 3)

        """
        x = ctx.convert(x)
        y, n = libmp.mpf_frexp(x._mpf_)
        return ctx.make_mpf(y), n

    def fneg(ctx, x, **kwargs):
        """
        Negates the number *x*, giving a floating-point result, optionally
        using a custom precision and rounding mode.

        See the documentation of :func:`fadd` for a detailed description
        of how to specify precision and rounding.

        **Examples**

        An mpmath number is returned::

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> fneg(2.5)
            mpf('-2.5')
            >>> fneg(-5+2j)
            mpc(real='5.0', imag='-2.0')

        Precise control over rounding is possible::

            >>> x = fadd(2, 1e-100, exact=True)
            >>> fneg(x)
            mpf('-2.0')
            >>> fneg(x, rounding='f')
            mpf('-2.0000000000000004')

        Negating with and without roundoff::

            >>> n = 200000000000000000000001
            >>> print int(-mpf(n))
            -200000000000000016777216
            >>> print int(fneg(n))
            -200000000000000016777216
            >>> print int(fneg(n, prec=log(n,2)+1))
            -200000000000000000000001
            >>> print int(fneg(n, dps=log(n,10)+1))
            -200000000000000000000001
            >>> print int(fneg(n, prec=inf))
            -200000000000000000000001
            >>> print int(fneg(n, dps=inf))
            -200000000000000000000001
            >>> print int(fneg(n, exact=True))
            -200000000000000000000001

        """
        prec, rounding = ctx._parse_prec(kwargs)
        x = ctx.convert(x)
        if hasattr(x, '_mpf_'):
            return ctx.make_mpf(mpf_neg(x._mpf_, prec, rounding))
        if hasattr(x, '_mpc_'):
            return ctx.make_mpc(mpc_neg(x._mpc_, prec, rounding))
        raise ValueError("Arguments need to be mpf or mpc compatible numbers")

    def fadd(ctx, x, y, **kwargs):
        """
        Adds the numbers *x* and *y*, giving a floating-point result,
        optionally using a custom precision and rounding mode.

        The default precision is the working precision of the context.
        You can specify a custom precision in bits by passing the *prec* keyword
        argument, or by providing an equivalent decimal precision with the *dps*
        keyword argument. If the precision is set to ``+inf``, or if the flag
        *exact=True* is passed, an exact addition with no rounding is performed.

        When the precision is finite, the optional *rounding* keyword argument
        specifies the direction of rounding. Valid options are ``'n'`` for
        nearest (default), ``'f'`` for floor, ``'c'`` for ceiling, ``'d'``
        for down, ``'u'`` for up.

        **Examples**

        Using :func:`fadd` with precision and rounding control::

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> fadd(2, 1e-20)
            mpf('2.0')
            >>> fadd(2, 1e-20, rounding='u')
            mpf('2.0000000000000004')
            >>> nprint(fadd(2, 1e-20, prec=100), 25)
            2.00000000000000000001
            >>> nprint(fadd(2, 1e-20, dps=15), 25)
            2.0
            >>> nprint(fadd(2, 1e-20, dps=25), 25)
            2.00000000000000000001
            >>> nprint(fadd(2, 1e-20, exact=True), 25)
            2.00000000000000000001

        Exact addition avoids cancellation errors, enforcing familiar laws
        of numbers such as `x+y-x = y`, which don't hold in floating-point
        arithmetic with finite precision::

            >>> x, y = mpf(2), mpf('1e-1000')
            >>> print x + y - x
            0.0
            >>> print fadd(x, y, prec=inf) - x
            1.0e-1000
            >>> print fadd(x, y, exact=True) - x
            1.0e-1000

        Exact addition can be inefficient and may be impossible to perform
        with large magnitude differences::

            >>> fadd(1, '1e-100000000000000000000', prec=inf)
            Traceback (most recent call last):
              ...
            OverflowError: the exact result does not fit in memory

        """
        prec, rounding = ctx._parse_prec(kwargs)
        x = ctx.convert(x)
        y = ctx.convert(y)
        try:
            if hasattr(x, '_mpf_'):
                if hasattr(y, '_mpf_'):
                    return ctx.make_mpf(mpf_add(x._mpf_, y._mpf_, prec, rounding))
                if hasattr(y, '_mpc_'):
                    return ctx.make_mpc(mpc_add_mpf(y._mpc_, x._mpf_, prec, rounding))
            if hasattr(x, '_mpc_'):
                if hasattr(y, '_mpf_'):
                    return ctx.make_mpc(mpc_add_mpf(x._mpc_, y._mpf_, prec, rounding))
                if hasattr(y, '_mpc_'):
                    return ctx.make_mpc(mpc_add(x._mpc_, y._mpc_, prec, rounding))
        except (ValueError, OverflowError):
            raise OverflowError(ctx._exact_overflow_msg)
        raise ValueError("Arguments need to be mpf or mpc compatible numbers")

    def fsub(ctx, x, y, **kwargs):
        """
        Subtracts the numbers *x* and *y*, giving a floating-point result,
        optionally using a custom precision and rounding mode.

        See the documentation of :func:`fadd` for a detailed description
        of how to specify precision and rounding.

        **Examples**

        Using :func:`fsub` with precision and rounding control::

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> fsub(2, 1e-20)
            mpf('2.0')
            >>> fsub(2, 1e-20, rounding='d')
            mpf('1.9999999999999998')
            >>> nprint(fsub(2, 1e-20, prec=100), 25)
            1.99999999999999999999
            >>> nprint(fsub(2, 1e-20, dps=15), 25)
            2.0
            >>> nprint(fsub(2, 1e-20, dps=25), 25)
            1.99999999999999999999
            >>> nprint(fsub(2, 1e-20, exact=True), 25)
            1.99999999999999999999

        Exact subtraction avoids cancellation errors, enforcing familiar laws
        of numbers such as `x-y+y = x`, which don't hold in floating-point
        arithmetic with finite precision::

            >>> x, y = mpf(2), mpf('1e1000')
            >>> print x - y + y
            0.0
            >>> print fsub(x, y, prec=inf) + y
            2.0
            >>> print fsub(x, y, exact=True) + y
            2.0

        Exact addition can be inefficient and may be impossible to perform
        with large magnitude differences::

            >>> fsub(1, '1e-100000000000000000000', prec=inf)
            Traceback (most recent call last):
              ...
            OverflowError: the exact result does not fit in memory

        """
        prec, rounding = ctx._parse_prec(kwargs)
        x = ctx.convert(x)
        y = ctx.convert(y)
        try:
            if hasattr(x, '_mpf_'):
                if hasattr(y, '_mpf_'):
                    return ctx.make_mpf(mpf_sub(x._mpf_, y._mpf_, prec, rounding))
                if hasattr(y, '_mpc_'):
                    return ctx.make_mpc(mpc_sub((x._mpf_, fzero), y._mpc_, prec, rounding))
            if hasattr(x, '_mpc_'):
                if hasattr(y, '_mpf_'):
                    return ctx.make_mpc(mpc_sub_mpf(x._mpc_, y._mpf_, prec, rounding))
                if hasattr(y, '_mpc_'):
                    return ctx.make_mpc(mpc_sub(x._mpc_, y._mpc_, prec, rounding))
        except (ValueError, OverflowError):
            raise OverflowError(ctx._exact_overflow_msg)
        raise ValueError("Arguments need to be mpf or mpc compatible numbers")

    def fmul(ctx, x, y, **kwargs):
        """
        Multiplies the numbers *x* and *y*, giving a floating-point result,
        optionally using a custom precision and rounding mode.

        See the documentation of :func:`fadd` for a detailed description
        of how to specify precision and rounding.

        **Examples**

        The result is an mpmath number::

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> fmul(2, 5.0)
            mpf('10.0')
            >>> fmul(0.5j, 0.5)
            mpc(real='0.0', imag='0.25')

        Avoiding roundoff::

            >>> x, y = 10**10+1, 10**15+1
            >>> print x*y
            10000000001000010000000001
            >>> print mpf(x) * mpf(y)
            1.0000000001e+25
            >>> print int(mpf(x) * mpf(y))
            10000000001000011026399232
            >>> print int(fmul(x, y))
            10000000001000011026399232
            >>> print int(fmul(x, y, dps=25))
            10000000001000010000000001
            >>> print int(fmul(x, y, exact=True))
            10000000001000010000000001

        Exact multiplication with complex numbers can be inefficient and may
        be impossible to perform with large magnitude differences between
        real and imaginary parts::

            >>> x = 1+2j
            >>> y = mpc(2, '1e-100000000000000000000')
            >>> fmul(x, y)
            mpc(real='2.0', imag='4.0')
            >>> fmul(x, y, rounding='u')
            mpc(real='2.0', imag='4.0000000000000009')
            >>> fmul(x, y, exact=True)
            Traceback (most recent call last):
              ...
            OverflowError: the exact result does not fit in memory

        """
        prec, rounding = ctx._parse_prec(kwargs)
        x = ctx.convert(x)
        y = ctx.convert(y)
        try:
            if hasattr(x, '_mpf_'):
                if hasattr(y, '_mpf_'):
                    return ctx.make_mpf(mpf_mul(x._mpf_, y._mpf_, prec, rounding))
                if hasattr(y, '_mpc_'):
                    return ctx.make_mpc(mpc_mul_mpf(y._mpc_, x._mpf_, prec, rounding))
            if hasattr(x, '_mpc_'):
                if hasattr(y, '_mpf_'):
                    return ctx.make_mpc(mpc_mul_mpf(x._mpc_, y._mpf_, prec, rounding))
                if hasattr(y, '_mpc_'):
                    return ctx.make_mpc(mpc_mul(x._mpc_, y._mpc_, prec, rounding))
        except (ValueError, OverflowError):
            raise OverflowError(ctx._exact_overflow_msg)
        raise ValueError("Arguments need to be mpf or mpc compatible numbers")

    def fdiv(ctx, x, y, **kwargs):
        """
        Divides the numbers *x* and *y*, giving a floating-point result,
        optionally using a custom precision and rounding mode.

        See the documentation of :func:`fadd` for a detailed description
        of how to specify precision and rounding.

        **Examples**

        The result is an mpmath number::

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> fdiv(3, 2)
            mpf('1.5')
            >>> fdiv(2, 3)
            mpf('0.66666666666666663')
            >>> fdiv(2+4j, 0.5)
            mpc(real='4.0', imag='8.0')

        The rounding direction and precision can be controlled::

            >>> fdiv(2, 3, dps=3)    # Should be accurate to at least 3 digits
            mpf('0.6666259765625')
            >>> fdiv(2, 3, rounding='d')
            mpf('0.66666666666666663')
            >>> fdiv(2, 3, prec=60)
            mpf('0.66666666666666667')
            >>> fdiv(2, 3, rounding='u')
            mpf('0.66666666666666674')

        Checking the error of a division by performing it at higher precision::

            >>> fdiv(2, 3) - fdiv(2, 3, prec=100)
            mpf('-3.7007434154172148e-17')

        Unlike :func:`fadd`, :func:`fmul`, etc., exact division is not
        allowed since the quotient of two floating-point numbers generally
        does not have an exact floating-point representation. (In the
        future this might be changed to allow the case where the division
        is actually exact.)

            >>> fdiv(2, 3, exact=True)
            Traceback (most recent call last):
              ...
            ValueError: division is not an exact operation

        """
        prec, rounding = ctx._parse_prec(kwargs)
        if not prec:
            raise ValueError("division is not an exact operation")
        x = ctx.convert(x)
        y = ctx.convert(y)
        if hasattr(x, '_mpf_'):
            if hasattr(y, '_mpf_'):
                return ctx.make_mpf(mpf_div(x._mpf_, y._mpf_, prec, rounding))
            if hasattr(y, '_mpc_'):
                return ctx.make_mpc(mpc_div((x._mpf_, fzero), y._mpc_, prec, rounding))
        if hasattr(x, '_mpc_'):
            if hasattr(y, '_mpf_'):
                return ctx.make_mpc(mpc_div_mpf(x._mpc_, y._mpf_, prec, rounding))
            if hasattr(y, '_mpc_'):
                return ctx.make_mpc(mpc_div(x._mpc_, y._mpc_, prec, rounding))
        raise ValueError("Arguments need to be mpf or mpc compatible numbers")

    def nint_distance(ctx, x):
        """
        Returns (n, d) where n is the nearest integer to x and d is the
        log-2 distance (i.e. distance in bits) of n from x. If d < 0,
        (-d) gives the bits of cancellation when n is subtracted from x.
        This function is intended to be used to check for cancellation
        at poles.
        """
        if hasattr(x, "_mpf_"):
            re = x._mpf_
            im_dist = ctx.ninf
        elif hasattr(x, "_mpc_"):
            re, im = x._mpc_
            isign, iman, iexp, ibc = im
            if iman:
                im_dist = iexp + ibc
            elif im == fzero:
                im_dist = ctx.ninf
            else:
                raise ValueError("requires a finite number")
        elif isinstance(x, int_types):
            return int(x), ctx.ninf
        elif isinstance(x, rational.mpq):
            p, q = x
            n, r = divmod(p, q)
            if 2*r >= q:
                n += 1
            elif not r:
                return n, ctx.ninf
            # log(p/q-n) = log((p-nq)/q) = log(p-nq) - log(q)
            d = bitcount(abs(p-n*q)) - bitcount(q)
            return n, d
        else:
            x = ctx.convert(x)
            if hasattr(x, "_mpf_") or hasattr(x, "_mpc_"):
                return ctx.nint_distance(x)
            else:
                raise TypeError("requires an mpf/mpc")
        sign, man, exp, bc = re
        shift = exp+bc
        if sign:
            man = -man
        if shift < -1:
            n = 0
            re_dist = shift
        elif man:
            if exp >= 0:
                n = man << exp
                re_dist = ctx.ninf
            else:
                if shift >= 0:
                    xfixed = man << shift
                else:
                    xfixed = man >> (-shift)
                n1 = xfixed >> bc
                n2 = -((-xfixed) >> bc)
                dist1 = abs(xfixed - (n1<<bc))
                dist2 = abs(xfixed - (n2<<bc))
                if dist1 < dist2:
                    re_dist = dist1
                    n = n1
                else:
                    re_dist = dist2
                    n = n2
                if re_dist:
                    re_dist = bitcount(re_dist) - bc
                else:
                    re_dist = ctx.ninf
        elif re == fzero:
            re_dist = ctx.ninf
            n = 0
        else:
            raise ValueError("requires a finite number")
        return n, max(re_dist, im_dist)

    def fprod(ctx, factors):
        r"""
        Calculates a product containing a finite number of factors (for
        infinite products, see :func:`nprod`). The factors will be
        converted to mpmath numbers.

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
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

            >>> from mpmath import *
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
                b = MPZ_ZERO
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

    def absmin(ctx, x):
        """
        Returns ``abs(x).a`` for an interval, or ``abs(x)`` for anything else.
        """
        if hasattr(x, '_mpi_'):
            return abs(x).a
        return abs(x)

    def absmax(ctx, x):
        """
        Returns ``abs(x).b`` for an interval, or ``abs(x)`` for anything else.
        """
        if hasattr(x, '_mpi_'):
            return abs(x).b
        return abs(x)

    def _as_points(ctx, x):
        if hasattr(x, '_mpi_'):
            return [x.a, x.b]
        return x

    '''
    def _zetasum(ctx, s, a, b):
        """
        Computes sum of k^(-s) for k = a, a+1, ..., b with a, b both small
        integers.
        """
        a = int(a)
        b = int(b)
        s = ctx.convert(s)
        prec, rounding = ctx._prec_rounding
        if hasattr(s, '_mpf_'):
            v = ctx.make_mpf(libmp.mpf_zetasum(s._mpf_, a, b, prec))
        elif hasattr(s, '_mpc_'):
            v = ctx.make_mpc(libmp.mpc_zetasum(s._mpc_, a, b, prec))
        return v
    '''

    def _zetasum_fast(ctx, s, a, n, derivatives=[0], reflect=False):
        if not (ctx.isint(a) and hasattr(s, "_mpc_")):
            raise NotImplementedError
        a = int(a)
        prec = ctx._prec
        xs, ys = libmp.mpc_zetasum(s._mpc_, a, n, derivatives, reflect, prec)
        xs = map(ctx.make_mpc, xs)
        ys = map(ctx.make_mpc, ys)
        return xs, ys


class PrecisionManager:
    def __init__(self, ctx, precfun, dpsfun, normalize_output=False):
        self.ctx = ctx
        self.precfun = precfun
        self.dpsfun = dpsfun
        self.normalize_output = normalize_output
    def __call__(self, f):
        def g(*args, **kwargs):
            orig = self.ctx.prec
            try:
                if self.precfun:
                    self.ctx.prec = self.precfun(self.ctx.prec)
                else:
                    self.ctx.dps = self.dpsfun(self.ctx.dps)
                if self.normalize_output:
                    v = f(*args, **kwargs)
                    if type(v) is tuple:
                        return tuple([+a for a in v])
                    return +v
                else:
                    return f(*args, **kwargs)
            finally:
                self.ctx.prec = orig
        g.__name__ = f.__name__
        g.__doc__ = f.__doc__
        return g
    def __enter__(self):
        self.origp = self.ctx.prec
        if self.precfun:
            self.ctx.prec = self.precfun(self.ctx.prec)
        else:
            self.ctx.dps = self.dpsfun(self.ctx.dps)
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.ctx.prec = self.origp
        return False


if __name__ == '__main__':
    import doctest
    doctest.testmod()
