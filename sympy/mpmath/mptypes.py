"""
This module defines the mpf, mpc classes, and standard functions for
operating with them.
"""
__docformat__ = 'plaintext'

from settings import (MP_BASE, MP_ONE, mp, prec_rounding, extraprec, extradps,
    workprec, workdps, int_types, repr_dps, round_floor, round_ceiling)

from libmpf import (\
    ComplexResult, to_pickable, from_pickable, normalize,
    from_int, from_float, from_str, to_int, to_float, to_str,
    from_rational, from_man_exp,
    fone, fzero, finf, fninf, fnan,
    mpf_abs, mpf_pos, mpf_neg, mpf_add, mpf_sub, mpf_mul, mpf_mul_int,
    mpf_div, mpf_rdiv_int, mpf_pow_int, mpf_mod,
    mpf_eq, mpf_cmp, mpf_lt, mpf_gt, mpf_le, mpf_ge,
    mpf_hash, mpf_rand
)

from libmpc import (\
    complex_to_str,
    mpc_abs, mpc_add, mpc_add_mpf, mpc_sub, mpc_sub_mpf, mpc_mul, mpc_mul_mpf,
    mpc_mul_int, mpc_div, mpc_div_mpf, mpc_pow, mpc_pow_mpf, mpc_pow_int
)

from libelefun import mpf_pow

from libmpi import (\
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

def convert_lossless(x, strings=True):
    """Attempt to convert x to an mpf or mpc losslessly. If x is an
    mpf or mpc, return it unchanged. If x is an int, create an mpf with
    sufficient precision to represent it exactly. If x is a str, just
    convert it to an mpf with the current working precision (perhaps
    this should be done differently...)"""
    if isinstance(x, mpnumeric): return x
    if isinstance(x, int_types): return make_mpf(from_int(x))
    if isinstance(x, float): return make_mpf(from_float(x))
    if isinstance(x, complex): return mpc(x)
    if strings and isinstance(x, basestring): return make_mpf(from_str(x, *prec_rounding))
    if hasattr(x, '_mpf_'): return make_mpf(x._mpf_)
    if hasattr(x, '_mpc_'): return make_mpc(x._mpc_)
    if hasattr(x, '_mpmath_'): return convert_lossless(x._mpmath_(*prec_rounding))
    raise TypeError("cannot create mpf from " + repr(x))

def try_convert_mpf_value(x, prec, rounding):
    if isinstance(x, float): return from_float(x)
    if hasattr(x, '_mpf_'): return x._mpf_
    if hasattr(x, '_mpmath_'):
        t = convert_lossless(x._mpmath_(prec, rounding))
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
        t = convert_lossless(x._mpmath_(prec, rounding))
        if isinstance(t, mpf):
            return t._mpf_
    raise TypeError("cannot create mpf from " + repr(x))

def mpf_convert_rhs(x):
    if isinstance(x, int_types): return from_int(x)
    if isinstance(x, float): return from_float(x)
    if isinstance(x, complex_types): return mpc(x)
    if hasattr(x, '_mpf_'): return x._mpf_
    if hasattr(x, '_mpmath_'):
        t = convert_lossless(x._mpmath_(*prec_rounding))
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
        return convert_lossless(x)
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
        other = convert_lossless(other, strings=False)
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

    def __repr__(self):
        return mpi_str(self._val, mp.prec)

    __str__ = __repr__

    def __eq__(self, other):
        if not isinstance(other, mpi):
            try:
                other = mpi(other)
            except:
                return NotImplemented
        return (self.a == other.a) and (self.b == other.b)

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
        if dps: prec = dps_to_prec(prec)
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
    """Return an mpf chosen randomly from [0, 1)."""
    return make_mpf(mpf_rand(mp.prec))

def isnan(x):
    if not isinstance(x, mpf):
        return False
    return x._mpf_ == fnan

def isinf(x):
    if not isinstance(x, mpf):
        return False
    return x._mpf_ in (finf, fninf)

def isint(x):
    if isinstance(x, int_types):
        return True
    try:
        x = convert_lossless(x)
    except:
        return False
    if isinstance(x, mpf):
        if isnan(x) or isinf(x):
            return False
        return x == int(x)
    return False

def absmin(x):
    """
    Returns abs(x).a for an interval, or abs(x) for anything else.
    """
    if isinstance(x, mpi):
        return abs(x).a
    return abs(x)

def absmax(x):
    """
    Returns abs(x).b for an interval, or abs(x) for anything else.
    """
    if isinstance(x, mpi):
        return abs(x).b
    return abs(x)

def AS_POINTS(x):
    if isinstance(x, mpi):
        return [x.a, x.b]
    return x

def fraction(p, q):
    """Given Python integers p, q, return a lazy mpf with value p/q.
    The value is updated with the precision.

        >>> mp.dps = 15
        >>> a = fraction(1,100)
        >>> b = mpf(1)/100
        >>> print a; print b
        0.01
        0.01
        >>> mp.dps = 30
        >>> print a; print b
        0.01
        0.0100000000000000002081668171172
        >>> mp.dps = 15
    """
    return constant(lambda prec, rnd: from_rational(p, q, prec, rnd),
        '%s/%s' % (p, q))

from operator import gt, lt

def arange(*args):
    """arange([a,] b[, dt]) -> list [a, a + dt, a + 2*dt, ..., x] with x < b"""
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
    linspace(a, b, n, [endpoint=False]) -> n evenly spaced samples from a to b

    Instead of a, b you can specify an mpi interval.
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
    """
    Determine whether the difference between s and t is smaller
    than a given epsilon.

    Both a maximum relative difference and a maximum difference
    ('epsilons') may be specified. The absolute difference is
    defined as |s-t| and the relative difference is defined
    as |s-t|/max(|s|, |t|).

    If only one epsilon is given, both are set to the same value.
    If none is given, both epsilons are set to 2**(-prec+m) where
    prec is the current working precision and m is a small integer.
    """
    t = convert_lossless(t)
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
    """Convert an mpf or mpc to a decimal string literal with n significant
    digits. The small default value for n is chosen to make this function
    useful for printing collections of numbers.

    If x is a list or tuple, the function is applied to each element.
    For unrecognized classes, this simply returns str(x).

    There is also a companion function nprint that prints the string
    instead of returning it.

        >>> nstr([+pi, ldexp(1,-500)])
        '[3.14159, 3.05494e-151]'
        >>> print([+pi, ldexp(1,-500)])
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
    """Print the result of nstr(x, n)."""
    print nstr(x, n)
