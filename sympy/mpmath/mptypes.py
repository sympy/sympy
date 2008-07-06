"""
This module defines the mpf, mpc classes, and standard functions for
operating with them.
"""
__docformat__ = 'plaintext'

__all__ = ["mpnumeric", "mpf", "mpc", "pi", "e", "ln2", "ln10",
  "j", "sqrt", "hypot", "exp", "log", "cos", "sin", "tan", "atan", "atan2",
  "power", "asin", "acos", "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
  "arg", "degree", "rand", "inf", "nan", "floor", "ceil", "isnan", "almosteq",
  "ldexp", "fraction", "nstr", "nprint", "mp", "extraprec",
  "extradps", "workprec", "workdps", "eps", "convert_lossless", "make_mpf",
  "make_mpc", "sec", "csc", "cot", "sech", "csch", "coth",
  "asec", "acsc", "acot", "asech", "acsch", "acoth", "arange",
  "ln", "log10", "frexp", "radians", "degrees", "modf"]

from lib import *
from libmpc import *

rounding_table = {
  'floor' : round_floor,
  'ceiling' : round_ceiling,
  'down' : round_down,
  'up' : round_up,
  'nearest' : round_nearest,
  'default' : round_nearest
}

reverse_rounding_table = {
  round_floor : 'floor',
  round_ceiling : 'ceiling',
  round_down : 'down',
  round_up : 'up',
  round_nearest : 'nearest'
}

class Context(object):

    __slots__ = ['trap_complex']

    def __repr__(self):
        lines = ["Mpmath settings:",
            ("  mp.prec = %s" % self.prec).ljust(30) + "[default: 53]",
            ("  mp.dps = %s" % self.dps).ljust(30) + "[default: 15]",
            ("  mp.rounding = '%s'" % self.rounding).ljust(30) + "[default: 'nearest']",
            ("  mp.trap_complex = %s" % self.trap_complex).ljust(30) + "[default: False]",
        ]
        return "\n".join(lines)

    def default(self):
        global gp, gd, gr
        gp = 53
        gd = 15
        gr = round_nearest
        self.trap_complex = False

    def set_prec(self, n):
        global gp, gd
        gp = max(1, int(n))
        gd = prec_to_dps(n)

    def set_dps(self, n):
        global gp, gd
        gp = dps_to_prec(n)
        gd = max(1, int(n))

    def set_rounding(self, s):
        global gr
        try:
            gr = rounding_table[s]
        except KeyError:
            raise ValueError(("unknown rounding mode: %s.\n" % s) +
                "Value must be one of: %r" % rounding_table.keys())

    prec = property(lambda self: gp, set_prec)
    dps = property(lambda self: gd, set_dps)
    rounding = property(lambda self: reverse_rounding_table[gr], set_rounding)

mp = Context()
mp.default()


class PrecisionManager:

    def __init__(self, precfun, dpsfun, normalize_output=False):
        self.precfun = precfun
        self.dpsfun = dpsfun
        self.normalize_output = normalize_output

    def __call__(self, f):
        def g(*args, **kwargs):
            orig = mp.prec
            try:
                if self.precfun:
                    mp.prec = self.precfun(mp.prec)
                else:
                    mp.dps = self.dpsfun(mp.dps)
                if self.normalize_output:
                    v = f(*args, **kwargs)
                    if type(v) is tuple:
                        return tuple([+a for a in v])
                    return +v
                else:
                    return f(*args, **kwargs)
            finally:
                mp.prec = orig
        g.__name__ = f.__name__
        g.__doc__ = f.__doc__
        return g

    def __enter__(self):
        self.origp = mp.prec
        if self.precfun:
            mp.prec = self.precfun(mp.prec)
        else:
            mp.dps = self.dpsfun(mp.dps)

    def __exit__(self, exc_type, exc_val, exc_tb):
        mp.prec = self.origp
        return False

def extraprec(n, normalize_output=False):
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
    return PrecisionManager(lambda p: p + n, None, normalize_output)

def extradps(n, normalize_output=False):
    """
    This function is analogous to extraprec (see documentation)
    but changes the decimal precision instead of the number of bits.
    """
    return PrecisionManager(None, lambda d: d + n, normalize_output)

def workprec(n, normalize_output=False):
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
    return PrecisionManager(lambda p: n, None, normalize_output)

def workdps(n, normalize_output=False):
    """
    This function is analogous to workprec (see documentation)
    but changes the decimal precision instead of the number of bits.
    """
    return PrecisionManager(None, lambda d: n, normalize_output)

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

def convert_lossless(x):
    """Attempt to convert x to an mpf or mpc losslessly. If x is an
    mpf or mpc, return it unchanged. If x is an int, create an mpf with
    sufficient precision to represent it exactly. If x is a str, just
    convert it to an mpf with the current working precision (perhaps
    this should be done differently...)"""
    if isinstance(x, mpnumeric): return x
    if isinstance(x, int_types): return make_mpf(from_int(x))
    if isinstance(x, float): return make_mpf(from_float(x))
    if isinstance(x, complex): return mpc(x)
    if isinstance(x, basestring): return make_mpf(from_str(x, gp, gr))
    if hasattr(x, '_mpf_'): return make_mpf(x._mpf_)
    if hasattr(x, '_mpc_'): return make_mpc(x._mpc_)
    raise TypeError("cannot create mpf from " + repr(x))

def mpf_convert_arg(x):
    if isinstance(x, int_types): return from_int(x)
    if isinstance(x, float): return from_float(x)
    if isinstance(x, basestring): return from_str(x, gp, gr)
    if hasattr(x, '_mpf_'): return x._mpf_
    raise TypeError("cannot create mpf from " + repr(x))

def mpf_convert_rhs(x):
    if isinstance(x, int_types): return from_int(x)
    if isinstance(x, float): return from_float(x)
    if hasattr(x, '_mpf_'): return x._mpf_
    return NotImplemented

def mpf_convert_lhs(x):
    if isinstance(x, complex_types): return mpc(x)
    if isinstance(x, int_types): return make_mpf(from_int(x))
    if isinstance(x, float): return make_mpf(from_float(x))
    if hasattr(x, '_mpf_'): return make_mpf(x._mpf_)
    return NotImplemented

new = object.__new__

class mpf(mpnumeric):
    """An mpf instance holds a real-valued floating-point number. mpf:s
    work analogously to Python floats, but support arbitrary-precision
    arithmetic."""

    __slots__ = ['_mpf_']

    def __new__(cls, val=fzero):
        """A new mpf can be created from a Python float, an int, a
        or a decimal string representing a number in floating-point
        format."""
        if type(val) is cls:
            sign, man, exp, bc = val._mpf_
            if (not man) and exp:
                return val
            return make_mpf(normalize(sign, man, exp, bc, gp, gr))
        elif type(val) is tuple:
            if len(val) == 2:
                return make_mpf(from_man_exp(val[0], val[1], gp, gr))
            if len(val) == 4:
                sign, man, exp, bc = val
                return make_mpf(normalize(sign, MP_BASE(man), exp, bc, gp, gr))
            raise ValueError
        else:
            return make_mpf(fpos(mpf_convert_arg(val), gp, gr))

    man_exp = property(lambda self: self._mpf_[1:3])
    man = property(lambda self: self._mpf_[1])
    exp = property(lambda self: self._mpf_[2])
    bc = property(lambda self: self._mpf_[3])

    real = property(lambda self: self)
    imag = property(lambda self: zero)

    def __getstate__(self): return to_pickable(self._mpf_)
    def __setstate__(self, val): self._mpf_ = from_pickable(val)

    def __repr__(s): return "mpf('%s')" % to_str(s._mpf_, repr_dps(gp))
    def __str__(s): return to_str(s._mpf_, gd)
    def __hash__(s): return fhash(s._mpf_)
    def __int__(s): return long(to_int(s._mpf_))
    __long__ = __int__
    def __float__(s): return to_float(s._mpf_)
    def __complex__(s): return float(s) + 0j

    def __nonzero__(s):
        return s._mpf_ != fzero

    def __abs__(s): return make_mpf(fabs(s._mpf_, gp, gr))
    def __pos__(s): return make_mpf(fpos(s._mpf_, gp, gr))
    def __neg__(s): return make_mpf(fneg(s._mpf_, gp, gr))

    def __eq__(s, t):
        if type(t) is mpf:
            return feq(s._mpf_, t._mpf_)
        u = mpf_convert_rhs(t)
        if u is NotImplemented:
            if isinstance(t, complex_types):
                return (not t.imag) and s == t.real
            return u
        return feq(s._mpf_, u)

    def __ne__(s, t):
        v = s.__eq__(t)
        if v is NotImplemented:
            return v
        return not v

    def __cmp__(s, t):
        if type(t) is mpf:
            t = t._mpf_
        else:
            t = mpf_convert_rhs(t)
            if t is NotImplemented:
                return t
        return fcmp(s._mpf_, t)

    def __lt__(s, t):
        if type(t) is mpf:
            t = t._mpf_
        else:
            t = mpf_convert_rhs(t)
            if t is NotImplemented:
                return t
        return flt(s._mpf_, t)

    def __gt__(s, t):
        if type(t) is mpf:
            t = t._mpf_
        else:
            t = mpf_convert_rhs(t)
            if t is NotImplemented:
                return t
        return fgt(s._mpf_, t)

    def __le__(s, t):
        if type(t) is mpf:
            t = t._mpf_
        else:
            t = mpf_convert_rhs(t)
            if t is NotImplemented:
                return t
        return fle(s._mpf_, t)

    def __ge__(s, t):
        if type(t) is mpf:
            t = t._mpf_
        else:
            t = mpf_convert_rhs(t)
            if t is NotImplemented:
                return t
        return fge(s._mpf_, t)

    def binop(s, t, f):
        if isinstance(t, complex_types):
            s = mpc(s)
            if f is feq: return s == t
            if f is fmul: return s * t
            if f is fadd: return s + t
            if f is fsub: return s - t
            if f is fdiv: return s / t
            raise ValueError("bad operation")
        t = mpf_convert_rhs(t)
        if t is NotImplemented:
            return t
        return make_mpf(f(s._mpf_, t, gp, gr))

    def __add__(s, t):
        r = new(mpf)
        sval = s._mpf_
        try:
            r._mpf_ = fadd(sval, t._mpf_, gp, gr)
            return r
        except:
            if isinstance(t, int_types):
                r._mpf_ = fadd(sval, from_int(t), gp, gr)
                return r
            if isinstance(t, mpc):
                return t + s
            return s.binop(t, fadd)

    def __sub__(s, t):
        r = new(mpf)
        sval = s._mpf_
        try:
            r._mpf_ = fsub(sval, t._mpf_, gp, gr)
            return r
        except:
            if isinstance(t, int_types):
                r._mpf_ = fadd(sval, from_int(-t), gp, gr)
                return r
            # TODO  if isinstance(t, mpc)
            return s.binop(t, fsub)

    def __mul__(s, t):
        r = new(mpf)
        sval = s._mpf_
        try:
            r._mpf_ = fmul(sval, t._mpf_, gp, gr)
            return r
        except:
            if isinstance(t, int_types):
                r._mpf_ = fmuli(sval, t, gp, gr)
                return r
            if isinstance(t, mpc):
                return t * s
            return s.binop(t, fmul)

    def __rmul__(s, t):
        if isinstance(t, int_types):
            r = new(mpf)
            r._mpf_ = fmuli(s._mpf_, t, gp, gr)
            return r
        return s.binop(t, fmul)

    def __div__(s, t):
        if isinstance(t, mpf):
            return make_mpf(fdiv(s._mpf_, t._mpf_, gp, gr))
        return s.binop(t, fdiv)

    def __mod__(s, t):
        if isinstance(t, mpf):
            return make_mpf(fmod(s._mpf_, t._mpf_, gp, gr))
        return s.binop(t, fmod)

    def __pow__(s, t):
        if isinstance(t, int_types):
            return make_mpf(fpowi(s._mpf_, t, gp, gr))
        if not isinstance(t, mpf):
            if isinstance(t, complex_types):
                t = convert_lossless(t)
                return make_mpc(mpc_pow((s._mpf_, fzero), t._mpc_, gp, gr))
            t = mpf_convert_rhs(t)
            if t is NotImplemented:
                return t
        else:
            t = t._mpf_

        try:
            return make_mpf(fpow(s._mpf_, t, gp, gr))
        except ComplexResult:
            if mp.trap_complex:
                raise
            return make_mpc(mpc_pow((s._mpf_, fzero), (t, fzero), gp, gr))

    __radd__ = __add__

    def __rsub__(s, t): return mpf_convert_lhs(t) - s

    def __rdiv__(s, t):
        if isinstance(t, int_types):
            return make_mpf(fdivi(t, s._mpf_, gp, gr))
        return mpf_convert_lhs(t) / s

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __rpow__(s, t): return mpf_convert_lhs(t) ** s
    def __rmod__(s, t): return mpf_convert_lhs(t) % s

    def sqrt(s):
        return sqrt(s)

    def ae(s, t, rel_eps=None, abs_eps=None):
        return almosteq(s, t, rel_eps, abs_eps)


def make_mpf(v, construct=object.__new__, cls=mpf):
    """Create mpf verbatim from a given tuple of data."""
    a = construct(cls)
    a._mpf_ = v
    return a

one = make_mpf(fone)
zero = make_mpf(fzero)
inf = make_mpf(finf)
ninf = make_mpf(fninf)
nan = make_mpf(fnan)

def isnan(x):
    if not isinstance(x, mpf):
        return False
    return x._mpf_ == fnan

class mpc(mpnumeric):
    """An mpc represents a complex number using a pair of mpf:s (one
    for the real part and another for the imaginary part.) The mpc
    class behaves fairly similarly to Python's complex type."""

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
        return "(%s)" % complex_to_str(s.real._mpf_, s.imag._mpf_, gd)

    def __complex__(s): return complex(float(s.real), float(s.imag))
    def __pos__(s): return mpc(s.real, s.imag)
    def __abs__(s): return make_mpf(mpc_abs(s._mpc_, gp, gr))
    def __neg__(s): return mpc(-s.real, -s.imag)
    def __nonzero__(s): return bool(s.real) or bool(s.imag)
    def conjugate(s): return mpc(s.real, -s.imag)

    def __eq__(s, t):
        if not isinstance(t, mpc):
            if isinstance(t, str):
                return False
            t = mpc(t)
        return s.real == t.real and s.imag == t.imag

    def __ne__(s, t): return not s.__eq__(t)

    def _compare(*args):
        raise TypeError("no ordering relation is defined for complex numbers")

    __gt__ = _compare
    __le__ = _compare
    __gt__ = _compare
    __ge__ = _compare

    def __add__(s, t):
        if not isinstance(t, mpc):
            try:
                t = mpf(t)
                return make_mpc(mpc_add_mpf(s._mpc_, t._mpf_, gp, gr))
            except:
                t = mpc(t)
        return make_mpc(mpc_add(s._mpc_, t._mpc_, gp, gr))

    def __sub__(s, t):
        if not isinstance(t, mpc):
            try:
                t = mpf(t)
                return make_mpc(mpc_sub_mpf(s._mpc_, t._mpf_, gp, gr))
            except:
                t = mpc(t)
        return make_mpc(mpc_sub(s._mpc_, t._mpc_, gp, gr))

    def __mul__(s, t):
        if not isinstance(t, mpc):
            if isinstance(t, int_types):
                return make_mpc(mpc_mul_int(s._mpc_, t, gp, gr))
            if isinstance(t, mpf):
                return make_mpc(mpc_mul_mpf(s._mpc_, t._mpf_, gp, gr))
            t = mpc(t)
        return make_mpc(mpc_mul(s._mpc_, t._mpc_, gp, gr))

    def __rmul__(s, t):
        if isinstance(t, int_types):
            return make_mpc(mpc_mul_int(s._mpc_, t, gp, gr))
        if not isinstance(t, mpc):
            t = mpc(t)
        return make_mpc(mpc_mul(s._mpc_, t._mpc_, gp, gr))

    def __div__(s, t):
        if not isinstance(t, mpc):
            try:
                t = mpf(t)
                return make_mpc(mpc_div_mpf(s._mpc_, t._mpf_, gp, gr))
            except:
                t = mpc(t)
        return make_mpc(mpc_div(s._mpc_, t._mpc_, gp, gr))

    def __pow__(s, t):
        if isinstance(t, int_types):
            return make_mpc(mpc_pow_int(s._mpc_, t, gp, gr))
        t = convert_lossless(t)
        if isinstance(t, mpf):
            return make_mpc(mpc_pow_mpf(s._mpc_, t._mpf_, gp, gr))
        return make_mpc(mpc_pow(s._mpc_, t._mpc_, gp, gr))

    __radd__ = __add__

    def __rsub__(s, t): return (-s) + t
    def __rpow__(s, t): return convert_lossless(t) ** s
    def __rdiv__(s, t): return mpc(t) / s

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def ae(s, t, rel_eps=None, abs_eps=None):
        return almosteq(s, t, rel_eps, abs_eps)


complex_types = (complex, mpc)

def make_mpc(v, construct=object.__new__, cls=mpc):
    a = construct(cls)
    a._mpc_ = v
    return a

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

    @property
    def _mpf_(self):
        return self.func(gp, gr)

    def __repr__(self):
        return "<%s: %s~>" % (self.name, nstr(self))

pi = constant(fpi, "pi")
degree = constant(fdegree, "degree")
e = constant(fe, "e")
ln2 = constant(flog2, "log 2")
ln10 = constant(flog10, "log 10")
eps = constant(lambda p, r: (0, MP_ONE, -p+1, 1), "epsilon of working precision")

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

def mpfunc(name, real_f, complex_f, doc):
    def f(x):
        prec, rnd = gp, gr
        if not isinstance(x, mpnumeric):
            x = convert_lossless(x)
        if isinstance(x, mpf):
            try:
                return make_mpf(real_f(x._mpf_, prec, rnd))
            except ComplexResult:
                if mp.trap_complex:
                    raise
                x = (x._mpf_, fzero)
        else:
            x = x._mpc_
        return make_mpc(complex_f(x, prec, rnd))
    f.__name__ = name
    f.__doc__ = "Returns the %s of x" % doc
    return f

def altfunc(f, name, desc):
    def g(x):
        orig = mp.prec
        try:
            mp.prec = orig + 10
            return one/f(x)
        finally:
            mp.prec = orig
    g.__name__ = name
    g.__doc__ = "Returns the %s of x, 1/%s(x)" % (desc, f.__name__)
    return g

def altinvfunc(f, name, desc):
    def g(x):
        orig = mp.prec
        try:
            mp.prec = orig + 10
            return f(one/x)
        finally:
            mp.prec = orig
    g.__name__ = name
    g.__doc__ = "Returns the inverse %s of x, %s(1/x)" % (desc, f.__name__)
    return g

sqrt = mpfunc('sqrt', fsqrt, mpc_sqrt, "principal square root")
exp = mpfunc('exp', fexp, mpc_exp, "exponential function")
ln = mpfunc('ln', flog, mpc_log, "natural logarithm")

cos = mpfunc('cos', fcos, mpc_cos, "cosine")
sin = mpfunc('sin', fsin, mpc_sin, "sine")
tan = mpfunc('tan', ftan, mpc_tan, "tangent")
cosh = mpfunc('cosh', fcosh, mpc_cosh, "hyperbolic cosine")
sinh = mpfunc('sinh', fsinh, mpc_sinh, "hyperbolic sine")
tanh = mpfunc('tanh', ftanh, mpc_tanh, "hyperbolic tangent")

acos = mpfunc('acos', facos, mpc_acos, "inverse cosine")
asin = mpfunc('asin', fasin, mpc_asin, "inverse sine")
atan = mpfunc('atan', fatan, mpc_atan, "inverse tangent")
asinh = mpfunc('asinh', fasinh, mpc_asinh, "inverse hyperbolic sine")
acosh = mpfunc('acosh', facosh, mpc_acosh, "inverse hyperbolic cosine")
atanh = mpfunc('atanh', fatanh, mpc_atanh, "inverse hyperbolic tangent")

sec = altfunc(cos, 'sec', 'secant')
csc = altfunc(sin, 'csc', 'cosecant')
cot = altfunc(tan, 'cot', 'cotangent')
sech = altfunc(cosh, 'sech', 'hyperbolic secant')
csch = altfunc(sinh, 'csch', 'hyperbolic cosecant')
coth = altfunc(tanh, 'coth', 'hyperbolic cotangent')

asec = altinvfunc(acos, 'asec', 'secant')
acsc = altinvfunc(asin, 'acsc', 'cosecant')
acot = altinvfunc(atan, 'acot', 'cotangent')
asech = altinvfunc(acosh, 'asech', 'hyperbolic secant')
acsch = altinvfunc(asinh, 'acsch', 'hyperbolic cosecant')
acoth = altinvfunc(atanh, 'acoth', 'hyperbolic cotangent')

def hypot(x, y):
    """Returns the Euclidean distance sqrt(x*x + y*y). Both x and y
    must be real."""
    x = convert_lossless(x)
    y = convert_lossless(y)
    return make_mpf(fhypot(x._mpf_, y._mpf_, gp, gr))

def floor(x):
    """Computes the floor function of x. Note: returns an mpf, not a
    Python int. If x is larger than the precision, it will be rounded,
    not necessarily in the floor direction."""
    x = convert_lossless(x)
    return make_mpf(ffloor(x._mpf_, gp, gr))

def ceil(x):
    """Computes the ceiling function of x. Note: returns an mpf, not a
    Python int. If x is larger than the precision, it will be rounded,
    not necessarily in the ceiling direction."""
    x = convert_lossless(x)
    return make_mpf(fceil(x._mpf_, gp, gr))

def ldexp(x, n):
    """Calculate mpf(x) * 2**n efficiently. No rounding is performed."""
    x = convert_lossless(x)
    return make_mpf(fshift(x._mpf_, n))

def frexp(x):
    """Convert x to a scaled number y in the range [0.5, 1). Returns
    (y, n) such that x = y * 2**n. No rounding is performed."""
    x = convert_lossless(x)
    y, n = mpf_frexp(x._mpf_)
    return make_mpf(y), n

@extraprec(5)
def arg(x):
    """Returns the complex argument (phase) of x. The returned value is
    an mpf instance. The argument is here defined to satisfy
    -pi < arg(x) <= pi. On the negative real half-axis, it is taken to
    be +pi."""
    x = mpc(x)
    return atan2(x.imag, x.real)

def log(x, b=None):
    """Returns the base-b logarithm of x. If b is unspecified, return
    the natural (base-e) logarithm. log(x, b) is defined as
    log(x)/log(b). log(0) raises ValueError.

    The natural logarithm is real if x > 0 and complex if x < 0 or if x
    is complex. The principal branch of the complex logarithm is chosen,
    for which Im(log(x)) = -pi < arg(x) <= pi. """
    if b is not None:
        mp.prec += 3
        a = ln(x) / ln(b)
        mp.prec -= 3
        return +a
    return ln(x)

def log10(x):
    """Base-10 logarithm. Equivalent to log(x,10)."""
    return log(x,10)

def power(x, y):
    """Converts x and y to mpf or mpc and returns x**y = exp(y*log(x))."""
    return convert_lossless(x) ** convert_lossless(y)

def modf(x,y):
    """Converts x and y to mpf or mpc and returns x % y"""
    x = convert_lossless(x)
    y = convert_lossless(y)
    return x % y

def degrees(x):
    """Convert x given in radians to degrees"""
    return x / degree

def radians(x):
    """Convert x given in degrees to radians"""
    return x * degree

def atan2(y,x):
    """atan2(y, x) has the same magnitude as atan(y/x) but accounts for
    the signs of y and x. (Defined for real x and y only.)"""
    x = convert_lossless(x)
    y = convert_lossless(y)
    return make_mpf(fatan2(y._mpf_, x._mpf_, gp, gr))

def rand():
    """Return an mpf chosen randomly from [0, 1)."""
    return make_mpf(frand(mp.prec))

from operator import gt, lt

def arange(*args):
    """arange([a,] b[, dt]) -> list [a, a + dt, a + 2*dt, ..., b]"""
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
        ##print i, t, op(t, b)
        if op(t, b):
            result.append(t)
        else:
            break
    return result

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
        rel_eps = abs_eps = make_mpf((0, 1, -gp+4, 1))
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
    return str(x)

def nprint(x, n=6):
    """Print the result of nstr(x, n)."""
    print nstr(x, n)
