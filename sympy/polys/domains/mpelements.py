"""Real and complex elements with built-in truncation. """

from sympy.polys.domains.domainelement import DomainElement

from sympy.mpmath.ctx_mp_python import PythonMPContext, _mpf, _mpc, _constant
from sympy.mpmath.libmp import (MPZ_ONE, fzero, fone, finf, fninf, fnan,
    round_nearest, mpf_mul, mpf_abs, mpf_lt, mpc_abs, repr_dps, int_types,
    from_int, from_float, from_str)

class RealElement(_mpf, DomainElement):
    """An element of a real domain. """

    __slots__ = ['__mpf__']

    def _set_mpf(self, val):
        prec, rounding = self.context._prec_rounding
        tol = self.context.tol

        if mpf_lt(mpf_abs(val, prec, rounding), tol):
            self.__mpf__ = fzero
        else:
            self.__mpf__ = val

    _mpf_ = property(lambda self: self.__mpf__, _set_mpf)

    def parent(self):
        return self.context._parent

class ComplexElement(_mpc, DomainElement):
    """An element of a complex domain. """

    __slots__ = ['__mpc__']

    def _set_mpc(self, val):
        prec, rounding = self.context._prec_rounding
        tol = self.context.tol

        # norm = mpc_abs(val, prec, rounding)
        # tol = mpf_max(tol, mpf_mul(norm, tol))

        re, im = val

        if mpf_lt(mpf_abs(re, prec, rounding), tol):
            re = fzero
        if mpf_lt(mpf_abs(im, prec, rounding), tol):
            im = fzero

        self.__mpc__ = (re, im)

    _mpc_ = property(lambda self: self.__mpc__, _set_mpc)

    def parent(self):
        return self.context._parent

new = object.__new__

class MPContext(PythonMPContext):

    def __init__(ctx, prec=53, dps=None, tol=None):
        ctx._prec_rounding = [prec, round_nearest]

        if dps is None:
            ctx._set_prec(prec)
        else:
            ctx._set_dps(dps)

        ctx.mpf = type('RealElement', (RealElement,), {})
        ctx.mpc = type('ComplexElement', (ComplexElement,), {})
        ctx.mpf._ctxdata = [ctx.mpf, new, ctx._prec_rounding]
        ctx.mpc._ctxdata = [ctx.mpc, new, ctx._prec_rounding]
        ctx.mpf.context = ctx
        ctx.mpc.context = ctx
        ctx.constant = type('constant', (_constant,), {})
        ctx.constant._ctxdata = [ctx.mpf, new, ctx._prec_rounding]
        ctx.constant.context = ctx

        ctx.types = [ctx.mpf, ctx.mpc, ctx.constant]
        ctx.trap_complex = True
        ctx.pretty = True

        if tol is None:
            hundred = (0, 25, 2, 5)
            eps = (0, MPZ_ONE, 1-ctx.prec, 1)
            ctx.tol = mpf_mul(hundred, eps)
        else:
            ctx.tol = ctx._convert_tol(tol)

        ctx.tolerance = ctx.make_mpf(ctx.tol)

        ctx.zero = ctx.make_mpf(fzero)
        ctx.one = ctx.make_mpf(fone)
        ctx.j = ctx.make_mpc((fzero, fone))
        ctx.inf = ctx.make_mpf(finf)
        ctx.ninf = ctx.make_mpf(fninf)
        ctx.nan = ctx.make_mpf(fnan)

    def _convert_tol(ctx, tol):
        if isinstance(tol, int_types):
            return from_int(tol)
        if isinstance(tol, float):
            return from_float(tol)
        prec, rounding = ctx._prec_rounding
        if isinstance(tol, basestring):
            return from_str(tol, prec, rounding)
        raise ValueError("expected a real number, got %s" % tol)

    @property
    def _repr_digits(ctx):
        return repr_dps(ctx._prec)

    @property
    def _str_digits(ctx):
        return ctx._dps
