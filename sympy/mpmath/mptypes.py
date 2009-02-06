"""
This module defines the mpf, mpc classes, and standard functions for
operating with them.
"""
__docformat__ = 'plaintext'

import re

from settings import (MP_BASE, MP_ONE, mp, prec_rounding, extraprec, extradps,
    workprec, workdps, int_types, repr_dps, round_floor, round_ceiling, dps_to_prec)

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
    complex_to_str,
    mpc_abs, mpc_add, mpc_add_mpf, mpc_sub, mpc_sub_mpf, mpc_mul, mpc_mul_mpf,
    mpc_mul_int, mpc_div, mpc_div_mpf, mpc_pow, mpc_pow_mpf, mpc_pow_int
)

from libelefun import mpf_pow

from libmpi import (
    mpi_mid, mpi_delta, mpi_str,
    mpi_abs, mpi_pos, mpi_neg, mpi_add, mpi_sub,
    mpi_mul, mpi_div, mpi_pow_int, mpi_pow
)

class mpnumeric(object):
    """Base class for mpf and mpc. Calling mpnumeric(x) returns an mpf
    if x can be converted to an mpf (if it is a float, int, mpf, ...),
    and an mpc if x is complex."""
    __slots__ = []
    def __new__(cls, val):
        # TODO: should maybe normalize here
        if isinstance(val, cls): return val
        if isinstance(val, complex): return mpc(val)
        return mpf(val)

get_complex = re.compile(r'^\(?(?P<re>[\+\-]?\d*\.?\d*(e[\+\-]?\d+)?)??'
                         r'(?P<im>[\+\-]?\d*\.?\d*(e[\+\-]?\d+)?j)?\)?$')
# TODO: add tests

def mpmathify(x, strings=True):
    """
    Converts *x* to an ``mpf`` or ``mpc``. If *x* is of type ``mpf``,
    ``mpc``, ``int``, ``float``, ``complex``, the conversion
    will be performed losslessly.

    If *x* is a string, the result will be rounded to the present
    working precision. Strings representing fractions or complex
    numbers are permitted.

        >>> from sympy.mpmath import *
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
    if isinstance(x, mpnumeric): return x
    if isinstance(x, int_types): return make_mpf(from_int(x))
    if isinstance(x, float): return make_mpf(from_float(x))
    if isinstance(x, complex): return mpc(x)
    if strings and isinstance(x, basestring):
        try:
            return make_mpf(from_str(x, *prec_rounding))
        except Exception, e:
            if '/' in x:
                fract = x.split('/')
                assert len(fract) == 2
                return mpmathify(fract[0]) / mpmathify(fract[1])
            if 'j' in x.lower():
                x = x.lower().replace(' ', '')
                match = get_complex.match(x)
                re = match.group('re')
                if not re:
                    re = 0
                im = match.group('im').rstrip('j')
                return mpc(mpmathify(re),
                           mpmathify(im))
            raise e
    if hasattr(x, '_mpf_'): return make_mpf(x._mpf_)
    if hasattr(x, '_mpc_'): return make_mpc(x._mpc_)
    if hasattr(x, '_mpmath_'):
        return mpmathify(x._mpmath_(*prec_rounding))
    raise TypeError("cannot create mpf from " + repr(x))

def try_convert_mpf_value(x, prec, rounding):
    if isinstance(x, float): return from_float(x)
    if hasattr(x, '_mpf_'): return x._mpf_
    if hasattr(x, '_mpmath_'):
        t = mpmathify(x._mpmath_(prec, rounding))
        if isinstance(t, mpf):
            return t._mpf_
    return NotImplemented

def mpf_convert_arg(x, prec, rounding):
    if isinstance(x, int_types): return from_int(x)
    if isinstance(x, float): return from_float(x)
    if isinstance(x, basestring): return from_str(x, prec, rounding)
    if isinstance(x, constant): return x.func(prec, rounding)
    if hasattr(x, '_mpf_'): return x._mpf_
    if hasattr(x, '_mpmath_'):
        t = mpmathify(x._mpmath_(prec, rounding))
        if isinstance(t, mpf):
            return t._mpf_
    raise TypeError("cannot create mpf from " + repr(x))

def mpf_convert_rhs(x):
    if isinstance(x, int_types): return from_int(x)
    if isinstance(x, float): return from_float(x)
    if isinstance(x, complex_types): return mpc(x)
    if hasattr(x, '_mpf_'): return x._mpf_
    if hasattr(x, '_mpmath_'):
        t = mpmathify(x._mpmath_(*prec_rounding))
        if isinstance(t, mpf):
            return t._mpf_
        return t
    return NotImplemented

def mpf_convert_lhs(x):
    x = mpf_convert_rhs(x)
    if type(x) is tuple:
        return make_mpf(x)
    return x

def mpc_convert_lhs(x):
    try:
        return mpmathify(x)
    except TypeError:
        return NotImplemented

new = object.__new__

class mpf(mpnumeric):
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
        prec, rounding = prec_rounding
        if kwargs:
            prec = kwargs.get('prec', prec)
            if 'dps' in kwargs:
                prec = dps_to_prec(kwargs['dps'])
            rounding = kwargs.get('rounding', rounding)
        if type(val) is cls:
            sign, man, exp, bc = val._mpf_
            if (not man) and exp:
                return val
            return make_mpf(normalize(sign, man, exp, bc, prec, rounding))
        elif type(val) is tuple:
            if len(val) == 2:
                return make_mpf(from_man_exp(val[0], val[1], prec, rounding))
            if len(val) == 4:
                sign, man, exp, bc = val
                return make_mpf(normalize(sign, MP_BASE(man), exp, bc, prec, rounding))
            raise ValueError
        else:
            return make_mpf(mpf_pos(mpf_convert_arg(val, prec, rounding), prec, rounding))

    man_exp = property(lambda self: self._mpf_[1:3])
    man = property(lambda self: self._mpf_[1])
    exp = property(lambda self: self._mpf_[2])
    bc = property(lambda self: self._mpf_[3])

    real = property(lambda self: self)
    imag = property(lambda self: zero)

    conjugate = lambda self: self

    def __getstate__(self): return to_pickable(self._mpf_)
    def __setstate__(self, val): self._mpf_ = from_pickable(val)

    def __repr__(s): return "mpf('%s')" % to_str(s._mpf_, repr_dps(mp.prec))
    def __str__(s): return to_str(s._mpf_, mp.dps)
    def __hash__(s): return mpf_hash(s._mpf_)
    def __int__(s): return int(to_int(s._mpf_))
    def __long__(s): return long(to_int(s._mpf_))
    def __float__(s): return to_float(s._mpf_)
    def __complex__(s): return complex(float(s))
    def __nonzero__(s): return s._mpf_ != fzero
    def __abs__(s): return make_mpf(mpf_abs(s._mpf_, *prec_rounding))
    def __pos__(s): return make_mpf(mpf_pos(s._mpf_, *prec_rounding))
    def __neg__(s): return make_mpf(mpf_neg(s._mpf_, *prec_rounding))

    def _cmp(s, t, func):
        if hasattr(t, '_mpf_'):
            t = t._mpf_
        else:
            t = mpf_convert_rhs(t)
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
        prec, rounding = prec_rounding
        if type(t) in int_types:
            return make_mpf(mpf_sub(from_int(t), s._mpf_, prec, rounding))
        t = mpf_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t - s

    def __rdiv__(s, t):
        prec, rounding = prec_rounding
        if isinstance(t, int_types):
            return make_mpf(mpf_rdiv_int(t, s._mpf_, prec, rounding))
        t = mpf_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t / s

    def __rpow__(s, t):
        t = mpf_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t ** s

    def __rmod__(s, t):
        t = mpf_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t % s

    def sqrt(s):
        from functions import sqrt
        return sqrt(s)

    def ae(s, t, rel_eps=None, abs_eps=None):
        return almosteq(s, t, rel_eps, abs_eps)


mpf_binary_op = """
def %NAME%(self, other):
    prec, rounding = prec_rounding
    sval = self._mpf_
    if hasattr(other, '_mpf_'):
    #try:
        tval = other._mpf_
        %WITH_MPF%
    #except AttributeError:
    #    pass
    ttype = type(other)
    if ttype in int_types:
        %WITH_INT%
    elif ttype is float:
        tval = from_float(other)
        %WITH_MPF%
    elif ttype is mpc:
        tval = other._mpc_
        %WITH_MPC%
    elif ttype is complex:
        tval = from_float(other.real), from_float(other.imag)
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
            if mp.trap_complex:
                raise
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

mpf.__eq__ = binary_op('__eq__',
    'return mpf_eq(sval, tval)',
    'return mpf_eq(sval, from_int(other))',
    'return (tval[1] == fzero) and mpf_eq(tval[0], sval)')

mpf.__add__ = binary_op('__add__',
    'val = mpf_add(sval, tval, prec, rounding)' + return_mpf,
    'val = mpf_add(sval, from_int(other), prec, rounding)' + return_mpf,
    'val = mpc_add_mpf(tval, sval, prec, rounding)' + return_mpc)

mpf.__sub__ = binary_op('__sub__',
    'val = mpf_sub(sval, tval, prec, rounding)' + return_mpf,
    'val = mpf_sub(sval, from_int(other), prec, rounding)' + return_mpf,
    'val = mpc_sub((sval, fzero), tval, prec, rounding)' + return_mpc)

mpf.__mul__ = binary_op('__mul__',
    'val = mpf_mul(sval, tval, prec, rounding)' + return_mpf,
    'val = mpf_mul_int(sval, other, prec, rounding)' + return_mpf,
    'val = mpc_mul_mpf(tval, sval, prec, rounding)' + return_mpc)

mpf.__div__ = binary_op('__div__',
    'val = mpf_div(sval, tval, prec, rounding)' + return_mpf,
    'val = mpf_div(sval, from_int(other), prec, rounding)' + return_mpf,
    'val = mpc_div((sval, fzero), tval, prec, rounding)' + return_mpc)

mpf.__mod__ = binary_op('__mod__',
    'val = mpf_mod(sval, tval, prec, rounding)' + return_mpf,
    'val = mpf_mod(sval, from_int(other), prec, rounding)' + return_mpf,
    'raise NotImplementedError("complex modulo")')

mpf.__pow__ = binary_op('__pow__',
    mpf_pow_same,
    'val = mpf_pow_int(sval, other, prec, rounding)' + return_mpf,
    'val = mpc_pow((sval, fzero), tval, prec, rounding)' + return_mpc)

mpf.__radd__ = mpf.__add__
mpf.__rmul__ = mpf.__mul__
mpf.__truediv__ = mpf.__div__
mpf.__rtruediv__ = mpf.__rdiv__


class mpc(mpnumeric):
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
        elif hasattr(real, "_mpc_"):
            s._mpc_ = real._mpc_
            return s
        real = mpf(real)
        imag = mpf(imag)
        s._mpc_ = (real._mpf_, imag._mpf_)
        return s

    real = property(lambda self: make_mpf(self._mpc_[0]))
    imag = property(lambda self: make_mpf(self._mpc_[1]))

    def __getstate__(self):
        return to_pickable(self._mpc_[0]), to_pickable(self._mpc_[1])

    def __setstate__(self, val):
        self._mpc_ = from_pickable(val[0]), from_pickable(val[1])

    def __repr__(s):
        r = repr(s.real)[4:-1]
        i = repr(s.imag)[4:-1]
        return "mpc(real=%s, imag=%s)" % (r, i)

    def __str__(s):
        return "(%s)" % complex_to_str(s.real._mpf_, s.imag._mpf_, mp.dps)

    def __complex__(s): return complex(float(s.real), float(s.imag))
    def __pos__(s): return mpc(s.real, s.imag)
    def __abs__(s): return make_mpf(mpc_abs(s._mpc_, *prec_rounding))
    def __neg__(s): return mpc(-s.real, -s.imag)
    def __nonzero__(s): return bool(s.real) or bool(s.imag)
    def conjugate(s): return mpc(s.real, -s.imag)

    def __eq__(s, t):
        if not isinstance(t, mpc):
            if isinstance(t, str):
                return False
            t = mpc_convert_lhs(t)
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
        prec, rounding = prec_rounding
        if not isinstance(t, mpc):
            t = mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
            if isinstance(t, mpf):
                return make_mpc(mpc_add_mpf(s._mpc_, t._mpf_, prec, rounding))
        return make_mpc(mpc_add(s._mpc_, t._mpc_, prec, rounding))

    def __sub__(s, t):
        prec, rounding = prec_rounding
        if not isinstance(t, mpc):
            t = mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
            if isinstance(t, mpf):
                return make_mpc(mpc_sub_mpf(s._mpc_, t._mpf_, prec, rounding))
        return make_mpc(mpc_sub(s._mpc_, t._mpc_, prec, rounding))

    def __mul__(s, t):
        prec, rounding = prec_rounding
        if not isinstance(t, mpc):
            if isinstance(t, int_types):
                return make_mpc(mpc_mul_int(s._mpc_, t, prec, rounding))
            t = mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
            if isinstance(t, mpf):
                return make_mpc(mpc_mul_mpf(s._mpc_, t._mpf_, prec, rounding))
            t = mpc(t)
        return make_mpc(mpc_mul(s._mpc_, t._mpc_, prec, rounding))

    def __div__(s, t):
        prec, rounding = prec_rounding
        if not isinstance(t, mpc):
            t = mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
            if isinstance(t, mpf):
                return make_mpc(mpc_div_mpf(s._mpc_, t._mpf_, prec, rounding))
        return make_mpc(mpc_div(s._mpc_, t._mpc_, prec, rounding))

    def __pow__(s, t):
        prec, rounding = prec_rounding
        if isinstance(t, int_types):
            return make_mpc(mpc_pow_int(s._mpc_, t, prec, rounding))
        t = mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        if isinstance(t, mpf):
            return make_mpc(mpc_pow_mpf(s._mpc_, t._mpf_, prec, rounding))
        return make_mpc(mpc_pow(s._mpc_, t._mpc_, prec, rounding))

    __radd__ = __add__

    def __rsub__(s, t):
        t = mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t - s

    def __rmul__(s, t):
        prec, rounding = prec_rounding
        if isinstance(t, int_types):
            return make_mpc(mpc_mul_int(s._mpc_, t, prec, rounding))
        t = mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t * s

    def __rdiv__(s, t):
        t = mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t / s

    def __rpow__(s, t):
        t = mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t ** s

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def ae(s, t, rel_eps=None, abs_eps=None):
        return almosteq(s, t, rel_eps, abs_eps)


complex_types = (complex, mpc)


class mpi(mpnumeric):
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
        return make_mpf(self._val[0])

    @property
    def b(self):
        return make_mpf(self._val[1])

    @property
    def mid(self):
        return make_mpf(mpi_mid(self._val, mp.prec))

    @property
    def delta(self):
        return make_mpf(mpi_delta(self._val, mp.prec))

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
        return mpi_str(self._val, mp.prec)

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
        return make_mpi(mpi_abs(self._val, mp.prec))

    def __pos__(self):
        return make_mpi(mpi_pos(self._val, mp.prec))

    def __neg__(self):
        return make_mpi(mpi_neg(self._val, mp.prec))

    def __add__(self, other):
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_add(self._val, other._val, mp.prec))

    def __sub__(self, other):
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_sub(self._val, other._val, mp.prec))

    def __mul__(self, other):
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_mul(self._val, other._val, mp.prec))

    def __div__(self, other):
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_div(self._val, other._val, mp.prec))

    def __pow__(self, other):
        if isinstance(other, (int, long)):
            return make_mpi(mpi_pow_int(self._val, int(other), mp.prec))
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_pow(self._val, other._val, mp.prec))

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

def make_mpi(val, cls=mpi):
    a = new(cls)
    a._val = val
    return a

def make_mpf(v, cls=mpf):
    a = new(cls)
    a._mpf_ = v
    return a

def make_mpc(v, cls=mpc):
    a = new(cls)
    a._mpc_ = v
    return a

one = make_mpf(fone)
zero = make_mpf(fzero)
inf = make_mpf(finf)
ninf = make_mpf(fninf)
nan = make_mpf(fnan)
j = mpc(0,1)

class constant(mpf):
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
        if not prec: prec = prec_rounding[0]
        if not rounding: rounding = prec_rounding[1]
        if dps: prec = dps_to_prec(dps)
        return make_mpf(self.func(prec, rounding))

    @property
    def _mpf_(self):
        prec, rounding = prec_rounding
        return self.func(prec, rounding)

    def __repr__(self):
        return "<%s: %s~>" % (self.name, nstr(self))

eps = constant(lambda prec, rnd: (0, MP_ONE, 1-prec, 1),
    "epsilon of working precision")


def rand():
    """
    Returns an ``mpf`` with value chosen randomly from `[0, 1)`.
    The number of randomly generated bits in the mantissa is equal
    to the working precision.
    """
    return make_mpf(mpf_rand(mp.prec))

def isnan(x):
    """
    For an ``mpf`` *x*, determines whether *x* is not-a-number (nan)::

        >>> from sympy.mpmath import *
        >>> isnan(nan), isnan(3)
        (True, False)
    """
    if not isinstance(x, mpf):
        return False
    return x._mpf_ == fnan

def isinf(x):
    """
    For an ``mpf`` *x*, determines whether *x* is infinite::

        >>> from sympy.mpmath import *
        >>> isinf(inf), isinf(-inf), isinf(3)
        (True, True, False)
    """
    if not isinstance(x, mpf):
        return False
    return x._mpf_ in (finf, fninf)

def isint(x):
    """
    For an ``mpf`` *x*, or any type that can be converted
    to ``mpf``, determines whether *x* is exactly
    integer-valued::

        >>> from sympy.mpmath import *
        >>> isint(3), isint(mpf(3)), isint(3.2)
        (True, True, False)
    """
    if isinstance(x, int_types):
        return True
    try:
        x = mpmathify(x)
    except:
        return False
    if isinstance(x, mpf):
        if isnan(x) or isinf(x):
            return False
        return x == int(x)
    return False

def absmin(x):
    """
    Returns ``abs(x).a`` for an interval, or ``abs(x)`` for anything else.
    """
    if isinstance(x, mpi):
        return abs(x).a
    return abs(x)

def absmax(x):
    """
    Returns ``abs(x).b`` for an interval, or ``abs(x)`` for anything else.
    """
    if isinstance(x, mpi):
        return abs(x).b
    return abs(x)

def AS_POINTS(x):
    if isinstance(x, mpi):
        return [x.a, x.b]
    return x

def fraction(p, q):
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
    return constant(lambda prec, rnd: from_rational(p, q, prec, rnd),
        '%s/%s' % (p, q))

from operator import gt, lt

def arange(*args):
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

        >>> from sympy.mpmath import *
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
    a, b, dt = mpf(a), mpf(b), mpf(dt)
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

def linspace(*args, **kwargs):
    """
    ``linspace(a, b, n)`` returns a list of `n` evenly spaced
    samples from `a` to `b`. The syntax ``linspace(mpi(a,b), n)``
    is also valid.

    This function is often more convenient than :func:`arange`
    for partitioning an interval into subintervals, since
    the endpoint is included::

        >>> from sympy.mpmath import *
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
        a = mpf(args[0])
        b = mpf(args[1])
        n = int(args[2])
    elif len(args) == 2:
        assert isinstance(args[0], mpi)
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
            return [mpf(a)]
        step = (b - a) / mpf(n - 1)
        y = [i*step + a for i in xrange(n)]
        y[-1] = b
    else:
        step = (b - a) / mpf(n)
        y = [i*step + a for i in xrange(n)]
    return y

def almosteq(s, t, rel_eps=None, abs_eps=None):
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

        >>> from sympy.mpmath import *
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
    t = mpmathify(t)
    if abs_eps is None and rel_eps is None:
        rel_eps = abs_eps = make_mpf((0, 1, -mp.prec+4, 1))
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

def nstr(x, n=6):
    """
    Convert an ``mpf`` or ``mpc`` to a decimal string literal with *n*
    significant digits. The small default value for *n* is chosen to
    make this function useful for printing collections of numbers
    (lists, matrices, etc).

    If *x* is a list or tuple, :func:`nstr` is applied recursively
    to each element. For unrecognized classes, :func:`nstr`
    simply returns ``str(x)``.

    The companion function :func:`nprint` prints the result
    instead of returning it.

        >>> from sympy.mpmath import *
        >>> nstr([+pi, ldexp(1,-500)])
        '[3.14159, 3.05494e-151]'
        >>> nprint([+pi, ldexp(1,-500)])
        [3.14159, 3.05494e-151]
    """
    if isinstance(x, list):
        return "[%s]" % (", ".join(nstr(c, n) for c in x))
    if isinstance(x, tuple):
        return "(%s)" % (", ".join(nstr(c, n) for c in x))
    if isinstance(x, mpf):
        return to_str(x._mpf_, n)
    if isinstance(x, mpc):
        return "(" + complex_to_str(x._mpc_[0], x._mpc_[1], n)  + ")"
    if isinstance(x, basestring):
        return repr(x)
    from matrices import matrix
    if isinstance(x, matrix):
        return x.__nstr__(n)
    return str(x)

def nprint(x, n=6):
    """
    Equivalent to ``print nstr(x, n)``.
    """
    print nstr(x, n)

def chop(x, tol=None):
    """
    Chops off small real or imaginary parts, or converts
    numbers close to zero to exact zeros. The input can be a
    single number or an iterable::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> chop(5+1e-10j, tol=1e-9)
        mpf('5.0')
        >>> nprint(chop([1.0, 1e-20, 3+1e-18j, -4, 2]))
        [1.0, 0.0, 3.0, -4.0, 2.0]

    The tolerance defaults to ``100*eps``.
    """
    if tol is None:
        tol = 100*eps
    try:
        x = mpmathify(x)
        absx = abs(x)
        if abs(x) < tol:
            return zero
        if isinstance(x, mpc):
            if abs(x.imag) < min(tol, absx*tol):
                return x.real
            if abs(x.real) < min(tol, absx*tol):
                return mpc(0, x.imag)
    except TypeError:
        if hasattr(x, "__iter__"):
            return [chop(a, tol) for a in x]
    return x

def monitor(f, input='print', output='print'):
    """
    Returns a wrapped copy of *f* that monitors evaluation by calling
    *input* with every input (*args*, *kwargs*) passed to *f* and
    *output* with every value returned from *f*. The default action
    (specify using the special string value ``'print'``) is to print
    inputs and outputs to stdout, along with the total evaluation
    count::

        >>> from sympy.mpmath import *
        >>> mp.dps = 5
        >>> diff(monitor(exp), 1)   # diff will eval f(x-h) and f(x+h)
        in  0 (mpf('0.99999999906867742538452148'),) {}
        out 0 mpf('2.7182818259274480055282064')
        in  1 (mpf('1.0000000009313225746154785'),) {}
        out 1 mpf('2.7182818309906424675501024')
        mpf('2.7182808')

    To disable either the input or the output handler, you may
    pass *None* as argument.

    Custom input and output handlers may be used e.g. to store
    results for later analysis::

        >>> mp.dps = 15
        >>> input = []
        >>> output = []
        >>> findroot(monitor(sin, input.append, output.append), 3.0)
        mpf('3.1415926535897932')
        >>> len(input)  # Count number of evaluations
        9
        >>> print input[3], output[3]
        ((mpf('3.1415076583334066'),), {}) 8.49952562843408e-5
        >>> print input[4], output[4]
        ((mpf('3.1415928201669122'),), {}) -1.66577118985331e-7

    """
    if not input:
        input = lambda v: None
    elif input == 'print':
        incount = [0]
        def input(value):
            args, kwargs = value
            print "in  %s %r %r" % (incount[0], args, kwargs)
            incount[0] += 1
    if not output:
        output = lambda v: None
    elif output == 'print':
        outcount = [0]
        def output(value):
            print "out %s %r" % (outcount[0], value)
            outcount[0] += 1
    def f_monitored(*args, **kwargs):
        input((args, kwargs))
        v = f(*args, **kwargs)
        output(v)
        return v
    return f_monitored

def fsum(terms):
    """
    Calculates a sum containing a finite number of terms (for infinite
    series, see :func:`nsum`). The terms will be converted to
    mpmath numbers. For len(terms) > 2, this function is generally
    faster and produces more accurate results than the builtin
    Python function :func:`sum`.

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> fsum([1, 2, 0.5, 7])
        mpf('10.5')
    """
    prec, rnd = prec_rounding
    real = []
    imag = []
    other = 0
    for term in terms:
        if hasattr(term, "_mpf_"):
            val = term._mpf_
        elif hasattr(term, "_mpc_"):
            val, im = term._mpc_
            imag.append(im)
        else:
            term = mpmathify(term)
            if hasattr(term, "_mpf_"):
                val = term._mpf_
            elif hasattr(term, "_mpc_"):
                val, im = term._mpc_
                imag.append(im)
            else:
                other += term
                continue
        real.append(val)
    s = mpf_sum(real, prec, rnd)
    if imag:
        s = make_mpc((s, mpf_sum(imag, prec, rnd)))
    else:
        s = make_mpf(s)
    if other is 0:
        return s
    else:
        return s + other

def fdot(A, B=None):
    r"""
    Computes the dot product of the iterables `A` and `B`,

    .. math ::

        \sum_{k=0} A_k B_k.

    Alternatively, :func:`fdot` accepts a single iterable of pairs.
    In other words, ``fdot(A,B)`` and ``fdot(zip(A,B))`` are equivalent.

    The elements are automatically converted to mpmath numbers.

    Examples::

        >>> from sympy.mpmath import *
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
    prec, rnd = prec_rounding
    real = []
    imag = []
    other = 0
    hasattr_ = hasattr
    types = (mpf, mpc)
    for a, b in A:
        if type(a) not in types: a = mpmathify(a)
        if type(b) not in types: b = mpmathify(b)
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
        s = make_mpc((s, mpf_sum(imag, prec, rnd)))
    else:
        s = make_mpf(s)
    if other is 0:
        return s
    else:
        return s + other

def fprod(factors):
    r"""
    Calculates a product containing a finite number of factors (for
    infinite products, see :func:`nprod`). The factors will be
    converted to mpmath numbers.

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> fprod([1, 2, 0.5, 7])
        mpf('7.0')

    """
    orig = mp.prec
    try:
        v = mpf(1)
        for p in factors:
            v *= p
    finally:
        mp.prec = orig
    return +v

def timing(f, *args, **kwargs):
    """
    Returns time elapsed for evaluating ``f()``. Optionally arguments
    may be passed to time the execution of ``f(*args, **kwargs)``.

    If the first call is very quick, ``f`` is called
    repeatedly and the best time is returned.
    """
    once = kwargs.get('once')
    if 'once' in kwargs:
        del kwargs['once']
    if args or kwargs:
        if len(args) == 1 and not kwargs:
            arg = args[0]
            g = lambda: f(arg)
        else:
            g = lambda: f(*args, **kwargs)
    else:
        g = f
    from timeit import default_timer as clock
    t1=clock(); v=g(); t2=clock(); t=t2-t1
    if t > 0.05 or once:
        return t
    for i in range(3):
        t1=clock();
        # Evaluate multiple times because the timer function
        # has a significant overhead
        g();g();g();g();g();g();g();g();g();g()
        t2=clock()
        t=min(t,(t2-t1)/10)
    return t

if __name__ == '__main__':
    import doctest
    doctest.testmod()
