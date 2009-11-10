"""
Adaptive numerical evaluation of SymPy expressions, using mpmath
for mathematical functions.
"""

from sympy.mpmath.libmpf import (from_int, from_rational, fzero, normalize,
        bitcount, round_nearest, to_str, fone, fnone, fhalf, from_float,
        to_float, fnone, to_int, mpf_lt, mpf_sqrt, mpf_cmp, mpf_abs,
        mpf_pow_int, mpf_shift, mpf_add, mpf_mul, mpf_neg)

import sympy.mpmath.libmpc as libmpc
from sympy.mpmath.settings import dps_to_prec
from sympy.mpmath import mpf, mpc, quadts, quadosc, mp, make_mpf
from sympy.mpmath.libelefun import mpf_pi, mpf_log, mpf_pow, mpf_sin, mpf_cos, \
        mpf_atan, mpf_atan2, mpf_e, mpf_exp
from sympy.mpmath.libmpf import MP_BASE, from_man_exp
from sympy.mpmath.calculus import nsum
from sympy.mpmath import inf as mpmath_inf

from sympy.mpmath.gammazeta import mpf_bernoulli

import math

from basic import Basic, C, S
from sympify import sympify

LG10 = math.log(10,2)

# Used in a few places as placeholder values to denote exponents and
# precision levels, e.g. of exact numbers. Must be careful to avoid
# passing these to mpmath functions or returning them in final results.
INF = 1e1000
MINUS_INF = -1e1000

# ~= 100 digits. Real men set this to INF.
DEFAULT_MAXPREC = 333

class PrecisionExhausted(ArithmeticError):
    pass

#----------------------------------------------------------------------------#
#                                                                            #
#              Helper functions for arithmetic and complex parts             #
#                                                                            #
#----------------------------------------------------------------------------#

"""
An mpf value tuple is a tuple of integers (sign, man, exp, bc)
representing a floating-point numbers.

A temporary result is a tuple (re, im, re_acc, im_acc) where
re and im are nonzero mpf value tuples representing approximate
numbers, or None to denote exact zeros.

re_acc, im_acc are integers denoting log2(e) where e is the estimated
relative accuracy of the respective complex part, but may be anything
if the corresponding complex part is None.

"""

def fastlog(x):
    """Fast approximation of log2(x) for an mpf value tuple x.

    Notes: Calculated as exponent + width of mantissa. This is an
    approximation for two reasons: 1) it gives the ceil(log2(abs(x)))
    value and 2) it is too high by 1 in the case that x is an exact
    power of 2. Although this is easy to remedy by testing to see if
    the odd mpf mantissa is 1 (indicating that one was dealing with
    an exact power of 2) that would decrease the speed and is not
    necessary as this is only being used as an approximation for the
    number of bits in x. The correct return value could be written as
    "x[2] + (x[3] if x[1]!=1 else 0)".

    """
    if not x or x == fzero:
        return MINUS_INF
    return x[2] + x[3]

def complex_accuracy(result):
    """
    Returns relative accuracy of a complex number with given accuracies
    for the real and imaginary parts. The relative accuracy is defined
    in the complex norm sense as ||z|+|error|| / |z| where error
    is equal to (real absolute error) + (imag absolute error)*i.

    The full expression for the (logarithmic) error can be approximated
    easily by using the max norm to approximate the complex norm.

    In the worst case (re and im equal), this is wrong by a factor
    sqrt(2), or by log2(sqrt(2)) = 0.5 bit.
    """
    re, im, re_acc, im_acc = result
    if not im:
        if not re:
            return INF
        return re_acc
    if not re:
        return im_acc
    re_size = fastlog(re)
    im_size = fastlog(im)
    absolute_error = max(re_size-re_acc, im_size-im_acc)
    relative_error = absolute_error - max(re_size, im_size)
    return -relative_error

def get_abs(expr, prec, options):
    re, im, re_acc, im_acc = evalf(expr, prec+2, options)
    if not re:
        re, re_acc, im, im_acc = im, im_acc, re, re_acc
    if im:
        return libmpc.mpc_abs((re, im), prec), None, re_acc, None
    else:
        return mpf_abs(re), None, re_acc, None

def get_complex_part(expr, no, prec, options):
    """no = 0 for real part, no = 1 for imaginary part"""
    workprec = prec
    i = 0
    while 1:
        res = evalf(expr, workprec, options)
        value, accuracy = res[no::2]
        if (not value) or accuracy >= prec:
            return value, None, accuracy, None
        workprec += max(30, 2**i)
        i += 1

def evalf_abs(expr, prec, options):
    return get_abs(expr.args[0], prec, options)

def evalf_re(expr, prec, options):
    return get_complex_part(expr.args[0], 0, prec, options)

def evalf_im(expr, prec, options):
    return get_complex_part(expr.args[0], 1, prec, options)

def finalize_complex(re, im, prec):
    assert re and im
    if re == fzero and im == fzero:
        raise ValueError("got complex zero with unknown accuracy")
    size_re = fastlog(re)
    size_im = fastlog(im)
    # Convert fzeros to scaled zeros
    if re == fzero:
        re = mpf_shift(fone, size_im-prec)
        size_re = fastlog(re)
    elif im == fzero:
        im = mpf_shift(fone, size_re-prec)
        size_im = fastlog(im)
    if size_re > size_im:
        re_acc = prec
        im_acc = prec + min(-(size_re - size_im), 0)
    else:
        im_acc = prec
        re_acc = prec + min(-(size_im - size_re), 0)
    return re, im, re_acc, im_acc

def chop_parts(value, prec):
    """
    Chop off tiny real or complex parts.
    """
    re, im, re_acc, im_acc = value
    # Method 1: chop based on absolute value
    if re and (fastlog(re) < -prec+4):
        re, re_acc = None, None
    if im and (fastlog(im) < -prec+4):
        im, im_acc = None, None
    # Method 2: chop if inaccurate and relatively small
    if re and im:
        delta = fastlog(re) - fastlog(im)
        if re_acc < 2 and (delta - re_acc <= -prec+4):
            re, re_acc = None, None
        if im_acc < 2 and (delta - im_acc >= prec-4):
            im, im_acc = None, None
    return re, im, re_acc, im_acc

def check_target(expr, result, prec):
    a = complex_accuracy(result)
    if a < prec:
        raise PrecisionExhausted("Failed to distinguish the expression: \n\n%s\n\n"
            "from zero. Try simplifying the input, using chop=True, or providing "
            "a higher maxprec for evalf" % (expr))

def get_integer_part(expr, no, options, return_ints=False):
    """
    With no = 1, computes ceiling(expr)
    With no = -1, computes floor(expr)

    Note: this function either gives the exact result or signals failure.
    """

    # The expression is likely less than 2^30 or so
    assumed_size = 30
    ire, iim, ire_acc, iim_acc = evalf(expr, assumed_size, options)

    # We now know the size, so we can calculate how much extra precision
    # (if any) is needed to get within the nearest integer
    if ire and iim:
        gap = max(fastlog(ire)-ire_acc, fastlog(iim)-iim_acc)
    elif ire:
        gap = fastlog(ire)-ire_acc
    elif iim:
        gap = fastlog(iim)-iim_acc
    else:
        # ... or maybe the expression was exactly zero
        return None, None, None, None

    margin = 10

    if gap >= -margin:
        ire, iim, ire_acc, iim_acc = evalf(expr, margin+assumed_size+gap, options)

    # We can now easily find the nearest integer, but to find floor/ceil, we
    # must also calculate whether the difference to the nearest integer is
    # positive or negative (which may fail if very close)
    def calc_part(expr, nexpr):
        nint = int(to_int(nexpr, round_nearest))
        expr = C.Add(expr, -nint, evaluate=False)
        x, _, x_acc, _ = evalf(expr, 10, options)
        check_target(expr, (x, None, x_acc, None), 3)
        nint += int(no*(mpf_cmp(x or fzero, fzero) == no))
        nint = from_int(nint)
        return nint, fastlog(nint) + 10

    re, im, re_acc, im_acc = None, None, None, None

    if ire:
        re, re_acc = calc_part(C.re(expr, evaluate=False), ire)
    if iim:
        im, im_acc = calc_part(C.im(expr, evaluate=False), iim)

    if return_ints:
        return int(to_int(re or fzero)), int(to_int(im or fzero))
    return re, im, re_acc, im_acc

def evalf_ceiling(expr, prec, options):
    return get_integer_part(expr.args[0], 1, options)

def evalf_floor(expr, prec, options):
    return get_integer_part(expr.args[0], -1, options)

#----------------------------------------------------------------------------#
#                                                                            #
#                            Arithmetic operations                           #
#                                                                            #
#----------------------------------------------------------------------------#

def add_terms(terms, prec, target_prec):
    """
    Helper for evalf_add. Adds a list of (mpfval, accuracy) terms.
    """
    if len(terms) == 1:
        if not terms[0]:
            # XXX: this is supposed to represent a scaled zero
            return mpf_shift(fone, target_prec), -1
        return terms[0]
    max_extra_prec = 2*prec
    sum_man, sum_exp, absolute_error = 0, 0, MINUS_INF
    for x, accuracy in terms:
        if not x:
            continue
        sign, man, exp, bc = x
        if sign:
            man = -man
        absolute_error = max(absolute_error, bc+exp-accuracy)
        delta = exp - sum_exp
        if exp >= sum_exp:
            # x much larger than existing sum?
            # first: quick test
            if (delta > max_extra_prec) and \
                ((not sum_man) or delta-bitcount(abs(sum_man)) > max_extra_prec):
                sum_man = man
                sum_exp = exp
            else:
                sum_man += (man << delta)
        else:
            delta = -delta
            # x much smaller than existing sum?
            if delta-bc > max_extra_prec:
                if not sum_man:
                    sum_man, sum_exp = man, exp
            else:
                sum_man = (sum_man << delta) + man
                sum_exp = exp
    if absolute_error == MINUS_INF:
        return None, None
    if not sum_man:
        # XXX: this is supposed to represent a scaled zero
        return mpf_shift(fone, absolute_error), -1
    if sum_man < 0:
        sum_sign = 1
        sum_man = -sum_man
    else:
        sum_sign = 0
    sum_bc = bitcount(sum_man)
    sum_accuracy = sum_exp + sum_bc - absolute_error
    r = normalize(sum_sign, sum_man, sum_exp, sum_bc, target_prec,
        round_nearest), sum_accuracy
    #print "returning", to_str(r[0],50), r[1]
    return r

def evalf_add(v, prec, options):
    args = v.args
    target_prec = prec
    i = 0

    oldmaxprec = options.get('maxprec', DEFAULT_MAXPREC)
    options['maxprec'] = min(oldmaxprec, 2*prec)

    try:
        while 1:
            terms = [evalf(arg, prec+10, options) for arg in args]
            re, re_acc = add_terms([(a[0],a[2]) for a in terms if a[0]], prec, target_prec)
            im, im_acc = add_terms([(a[1],a[3]) for a in terms if a[1]], prec, target_prec)
            accuracy = complex_accuracy((re, im, re_acc, im_acc))
            if accuracy >= target_prec:
                if options.get('verbose'):
                    print "ADD: wanted", target_prec, "accurate bits, got", re_acc, im_acc
                return re, im, re_acc, im_acc
            else:
                diff = target_prec - accuracy
                if (prec-target_prec) > options.get('maxprec', DEFAULT_MAXPREC):
                    return re, im, re_acc, im_acc

                prec = prec + max(10+2**i, diff)
                options['maxprec'] = min(oldmaxprec, 2*prec)
                if options.get('verbose'):
                    print "ADD: restarting with prec", prec
            i += 1
    finally:
        options['maxprec'] = oldmaxprec

def evalf_mul(v, prec, options):
    args = v.args
    # With guard digits, multiplication in the real case does not destroy
    # accuracy. This is also true in the complex case when considering the
    # total accuracy; however accuracy for the real or imaginary parts
    # separately may be lower.
    acc = prec
    target_prec = prec
    # XXX: big overestimate
    prec = prec + len(args) + 5
    direction = 0
    # Empty product is 1
    man, exp, bc = MP_BASE(1), 0, 1
    direction = 0
    complex_factors = []
    # First, we multiply all pure real or pure imaginary numbers.
    # direction tells us that the result should be multiplied by
    # i**direction
    for arg in args:
        re, im, re_acc, im_acc = evalf(arg, prec, options)
        if re and im:
            complex_factors.append((re, im, re_acc, im_acc))
            continue
        elif re:
            (s, m, e, b), w_acc = re, re_acc
        elif im:
            (s, m, e, b), w_acc = im, im_acc
            direction += 1
        else:
            return None, None, None, None
        direction += 2*s
        man *= m
        exp += e
        bc += b
        if bc > 3*prec:
            man >>= prec
            exp += prec
        acc = min(acc, w_acc)
    sign = (direction & 2) >> 1
    v = normalize(sign, man, exp, bitcount(man), prec, round_nearest)
    if complex_factors:
        # make existing real scalar look like an imaginary and
        # multiply by the remaining complex numbers
        re, im = v, (0, MP_BASE(0), 0, 0)
        for wre, wim, wre_acc, wim_acc in complex_factors:
            # acc is the overall accuracy of the product; we aren't
            # computing exact accuracies of the product.
            acc = min(acc,
                      complex_accuracy((wre, wim, wre_acc, wim_acc)))
            A = mpf_mul(re, wre, prec)
            B = mpf_mul(mpf_neg(im), wim, prec)
            C = mpf_mul(re, wim, prec)
            D = mpf_mul(im, wre, prec)
            re, xre_acc = add_terms([(A, acc), (B, acc)], prec, target_prec)
            im, xim_acc = add_terms([(C, acc), (D, acc)], prec, target_prec)

        if options.get('verbose'):
            print "MUL: wanted", target_prec, "accurate bits, got", acc
        # multiply by i
        if direction & 1:
            return mpf_neg(im), re, acc, acc
        else:
            return re, im, acc, acc
    else:
        # multiply by i
        if direction & 1:
            return None, v, None, acc
        else:
            return v, None, acc, None

def evalf_pow(v, prec, options):

    target_prec = prec
    base, exp = v.args

    # We handle x**n separately. This has two purposes: 1) it is much
    # faster, because we avoid calling evalf on the exponent, and 2) it
    # allows better handling of real/imaginary parts that are exactly zero
    if exp.is_Integer:
        p = exp.p
        # Exact
        if not p:
            return fone, None, prec, None
        # Exponentiation by p magnifies relative error by |p|, so the
        # base must be evaluated with increased precision if p is large
        prec += int(math.log(abs(p),2))
        re, im, re_acc, im_acc = evalf(base, prec+5, options)
        # Real to integer power
        if re and not im:
            return mpf_pow_int(re, p, target_prec), None, target_prec, None
        # (x*I)**n = I**n * x**n
        if im and not re:
            z = mpf_pow_int(im, p, target_prec)
            case = p % 4
            if case == 0: return z, None, target_prec, None
            if case == 1: return None, z, None, target_prec
            if case == 2: return mpf_neg(z), None, target_prec, None
            if case == 3: return None, mpf_neg(z), None, target_prec
        # Zero raised to an integer power
        if not re:
            return None, None, None, None
        # General complex number to arbitrary integer power
        re, im = libmpc.mpc_pow_int((re, im), p, prec)
        # Assumes full accuracy in input
        return finalize_complex(re, im, target_prec)

    # Pure square root
    if exp is S.Half:
        xre, xim, xre_acc, yim_acc = evalf(base, prec+5, options)
        # General complex square root
        if xim:
            re, im = libmpc.mpc_sqrt((xre or fzero, xim), prec)
            return finalize_complex(re, im, prec)
        if not xre:
            return None, None, None, None
        # Square root of a negative real number
        if mpf_lt(xre, fzero):
            return None, mpf_sqrt(mpf_neg(xre), prec), None, prec
        # Positive square root
        return mpf_sqrt(xre, prec), None, prec, None

    # We first evaluate the exponent to find its magnitude
    # This determines the working precision that must be used
    prec += 10
    yre, yim, yre_acc, yim_acc = evalf(exp, prec, options)
    # Special cases: x**0
    if not (yre or yim):
        return fone, None, prec, None

    ysize = fastlog(yre)
    # Restart if too big
    # XXX: prec + ysize might exceed maxprec
    if ysize > 5:
        prec += ysize
        yre, yim, yre_acc, yim_acc = evalf(exp, prec, options)

    # Pure exponential function; no need to evalf the base
    if base is S.Exp1:
        if yim:
            re, im = libmpc.mpc_exp((yre or fzero, yim), prec)
            return finalize_complex(re, im, target_prec)
        return mpf_exp(yre, target_prec), None, target_prec, None

    xre, xim, xre_acc, yim_acc = evalf(base, prec+5, options)
    # 0**y
    if not (xre or xim):
        return None, None, None, None

    # (real ** complex) or (complex ** complex)
    if yim:
        re, im = libmpc.mpc_pow((xre or fzero, xim or fzero), (yre or fzero, yim),
            target_prec)
        return finalize_complex(re, im, target_prec)
    # complex ** real
    if xim:
        re, im = libmpc.mpc_pow_mpf((xre or fzero, xim), yre, target_prec)
        return finalize_complex(re, im, target_prec)
    # negative ** real
    elif mpf_lt(xre, fzero):
        re, im = libmpc.mpc_pow_mpf((xre, fzero), yre, target_prec)
        return finalize_complex(re, im, target_prec)
    # positive ** real
    else:
        return mpf_pow(xre, yre, target_prec), None, target_prec, None




#----------------------------------------------------------------------------#
#                                                                            #
#                            Special functions                               #
#                                                                            #
#----------------------------------------------------------------------------#

def evalf_trig(v, prec, options):
    """
    This function handles sin and cos of real arguments.

    TODO: should also handle tan and complex arguments.
    """
    if v.func is C.cos:
        func = mpf_cos
    elif v.func is C.sin:
        func = mpf_sin
    else:
        raise NotImplementedError
    arg = v.args[0]
    # 20 extra bits is possibly overkill. It does make the need
    # to restart very unlikely
    xprec = prec + 20
    re, im, re_acc, im_acc = evalf(arg, xprec, options)
    if im:
        raise NotImplementedError
    if not re:
        if v.func is C.cos:
            return fone, None, prec, None
        elif v.func is C.sin:
            return None, None, None, None
        else:
            raise NotImplementedError
    # For trigonometric functions, we are interested in the
    # fixed-point (absolute) accuracy of the argument.
    xsize = fastlog(re)
    # Magnitude <= 1.0. OK to compute directly, because there is no
    # danger of hitting the first root of cos (with sin, magnitude
    # <= 2.0 would actually be ok)
    if xsize < 1:
        return func(re, prec, round_nearest), None, prec, None
    # Very large
    if xsize >= 10:
        xprec = prec + xsize
        re, im, re_acc, im_acc = evalf(arg, xprec, options)
    # Need to repeat in case the argument is very close to a
    # multiple of pi (or pi/2), hitting close to a root
    while 1:
        y = func(re, prec, round_nearest)
        ysize = fastlog(y)
        gap = -ysize
        accuracy = (xprec - xsize) - gap
        if accuracy < prec:
            if options.get('verbose'):
                print "SIN/COS", accuracy, "wanted", prec, "gap", gap
                print to_str(y,10)
            if xprec > options.get('maxprec', DEFAULT_MAXPREC):
                return y, None, accuracy, None
            xprec += gap
            re, im, re_acc, im_acc = evalf(arg, xprec, options)
            continue
        else:
            return y, None, prec, None

def evalf_log(expr, prec, options):
    arg = expr.args[0]
    workprec = prec+10
    xre, xim, xacc, _ = evalf(arg, workprec, options)

    if xim:
        # XXX: use get_abs etc instead
        re = evalf_log(C.log(C.abs(arg, evaluate=False), evaluate=False), prec, options)
        im = mpf_atan2(xim, xre or fzero, prec)
        return re[0], im, re[2], prec

    imaginary_term = (mpf_cmp(xre, fzero) < 0)

    re = mpf_log(mpf_abs(xre), prec, round_nearest)
    size = fastlog(re)
    if prec - size > workprec:
        # We actually need to compute 1+x accurately, not x
        arg = C.Add(S.NegativeOne,arg,evaluate=False)
        xre, xim, xre_acc, xim_acc = evalf_add(arg, prec, options)
        prec2 = workprec - fastlog(xre)
        re = mpf_log(mpf_add(xre, fone, prec2), prec, round_nearest)

    re_acc = prec

    if imaginary_term:
        return re, mpf_pi(prec), re_acc, prec
    else:
        return re, None, re_acc, None

def evalf_atan(v, prec, options):
    arg = v.args[0]
    xre, xim, reacc, imacc = evalf(arg, prec+5, options)
    if xim:
        raise NotImplementedError
    return mpf_atan(xre, prec, round_nearest), None, prec, None

def evalf_piecewise(expr, prec, options):
    if 'subs' in options:
        expr = expr.subs(options['subs'])
        del options['subs']
        if hasattr(expr,'func'):
            return evalf(expr, prec, options)
        if type(expr) == float:
            return evalf(C.Real(expr), prec, options)
        if type(expr) == int:
            return evalf(C.Integer(expr), prec, options)

    # We still have undefined symbols
    raise NotImplementedError

def evalf_piecewise(expr, prec, options):
    if 'subs' in options:
        expr = expr.subs(options['subs'])
        del options['subs']
        if hasattr(expr,'func'):
            return evalf(expr, prec, options)
        if type(expr) == float:
            return evalf(C.Real(expr), prec, options)
        if type(expr) == int:
            return evalf(C.Integer(expr), prec, options)

    # We still have undefined symbols
    raise NotImplementedError

def evalf_bernoulli(expr, prec, options):
    arg = expr.args[0]
    if not arg.is_Integer:
        raise ValueError("Bernoulli number index must be an integer")
    n = int(arg)
    b = mpf_bernoulli(n, prec, round_nearest)
    if b == fzero:
        return None, None, None, None
    return b, None, prec, None

#----------------------------------------------------------------------------#
#                                                                            #
#                            High-level operations                           #
#                                                                            #
#----------------------------------------------------------------------------#

def as_mpmath(x, prec, options):
    x = sympify(x)
    if isinstance(x, C.Zero):
        return mpf(0)
    if isinstance(x, C.Infinity):
        return mpf('inf')
    if isinstance(x, C.NegativeInfinity):
        return mpf('-inf')
    # XXX
    re, im, _, _ = evalf(x, prec, options)
    if im:
        return mpc(re or fzero, im)
    return mpf(re)

def do_integral(expr, prec, options):
    func = expr.args[0]
    x, (xlow, xhigh) = expr.args[1][0]
    orig = mp.prec

    oldmaxprec = options.get('maxprec', DEFAULT_MAXPREC)
    options['maxprec'] = min(oldmaxprec, 2*prec)

    try:
        mp.prec = prec+5
        xlow = as_mpmath(xlow, prec+15, options)
        xhigh = as_mpmath(xhigh, prec+15, options)

        # Integration is like summation, and we can phone home from
        # the integrand function to update accuracy summation style
        # Note that this accuracy is inaccurate, since it fails
        # to account for the variable quadrature weights,
        # but it is better than nothing

        have_part = [False, False]
        max_real_term = [MINUS_INF]
        max_imag_term = [MINUS_INF]

        def f(t):
            re, im, re_acc, im_acc = evalf(func, mp.prec, {'subs':{x:t}})

            have_part[0] = re or have_part[0]
            have_part[1] = im or have_part[1]

            max_real_term[0] = max(max_real_term[0], fastlog(re))
            max_imag_term[0] = max(max_imag_term[0], fastlog(im))

            if im:
                return mpc(re or fzero, im)
            return mpf(re or fzero)

        if options.get('quad') == 'osc':
            A = C.Wild('A', exclude=[x])
            B = C.Wild('B', exclude=[x])
            D = C.Wild('D')
            m = func.match(C.cos(A*x+B)*D)
            if not m:
                m = func.match(C.sin(A*x+B)*D)
            if not m:
                raise ValueError("An integrand of the form sin(A*x+B)*f(x) "
                  "or cos(A*x+B)*f(x) is required for oscillatory quadrature")
            period = as_mpmath(2*S.Pi/m[A], prec+15, options)
            result = quadosc(f, [xlow, xhigh], period=period)
            # XXX: quadosc does not do error detection yet
            quadrature_error = MINUS_INF
        else:
            result, quadrature_error = quadts(f, [xlow, xhigh], error=1)
            quadrature_error = fastlog(quadrature_error._mpf_)

    finally:
        options['maxprec'] = oldmaxprec
        mp.prec = orig

    if have_part[0]:
        re = result.real._mpf_
        if re == fzero:
            re = mpf_shift(fone, min(-prec,-max_real_term[0],-quadrature_error))
            re_acc = -1
        else:
            re_acc = -max(max_real_term[0]-fastlog(re)-prec, quadrature_error)
    else:
        re, re_acc = None, None

    if have_part[1]:
        im = result.imag._mpf_
        if im == fzero:
            im = mpf_shift(fone, min(-prec,-max_imag_term[0],-quadrature_error))
            im_acc = -1
        else:
            im_acc = -max(max_imag_term[0]-fastlog(im)-prec, quadrature_error)
    else:
        im, im_acc = None, None

    result = re, im, re_acc, im_acc
    return result

def evalf_integral(expr, prec, options):
    workprec = prec
    i = 0
    maxprec = options.get('maxprec', INF)
    while 1:
        result = do_integral(expr, workprec, options)
        accuracy = complex_accuracy(result)
        if accuracy >= prec or workprec >= maxprec:
            return result
        workprec += prec - max(-2**i, accuracy)
        i += 1

def check_convergence(numer, denom, n):
    """
    Returns (h, g, p) where
    -- h is:
        > 0 for convergence of rate 1/factorial(n)**h
        < 0 for divergence of rate factorial(n)**(-h)
        = 0 for geometric or polynomial convergence or divergence

    -- abs(g) is:
        > 1 for geometric convergence of rate 1/h**n
        < 1 for geometric divergence of rate h**n
        = 1 for polynomial convergence or divergence

        (g < 0 indicates an alternating series)

    -- p is:
        > 1 for polynomial convergence of rate 1/n**h
        <= 1 for polynomial divergence of rate n**(-h)

    """
    npol = C.Poly(numer, n)
    dpol = C.Poly(denom, n)
    p = npol.degree
    q = dpol.degree
    rate = q - p
    if rate:
        return rate, None, None
    constant = dpol.lead_term[0] / npol.lead_term[0]
    if abs(constant) != 1:
        return rate, constant, None
    if npol.degree == dpol.degree == 0:
        return rate, constant, 0
    pc = list(npol.iter_all_terms())[1][0]
    qc = list(dpol.iter_all_terms())[1][0]
    return rate, constant, qc-pc

def hypsum(expr, n, start, prec):
    """
    Sum a rapidly convergent infinite hypergeometric series with
    given general term, e.g. e = hypsum(1/factorial(n), n). The
    quotient between successive terms must be a quotient of integer
    polynomials.
    """
    from sympy import hypersimp, lambdify

    if start:
        expr = expr.subs(n, n+start)
    hs = hypersimp(expr, n)
    if hs is None:
        raise NotImplementedError("a hypergeometric series is required")
    num, den = hs.as_numer_denom()

    func1 = lambdify(n, num)
    func2 = lambdify(n, den)

    h, g, p = check_convergence(num, den, n)

    if h < 0:
        raise ValueError("Sum diverges like (n!)^%i" % (-h))

    # Direct summation if geometric or faster
    if h > 0 or (h == 0 and abs(g) > 1):
        one = MP_BASE(1) << prec
        term = expr.subs(n, 0)
        term = (MP_BASE(term.p) << prec) // term.q
        s = term
        k = 1
        while abs(term) > 5:
            term *= MP_BASE(func1(k-1))
            term //= MP_BASE(func2(k-1))
            s += term
            k += 1
        return from_man_exp(s, -prec)
    else:
        alt = g < 0
        if abs(g) < 1:
            raise ValueError("Sum diverges like (%i)^n" % abs(1/g))
        if p < 1 or (p == 1 and not alt):
            raise ValueError("Sum diverges like n^%i" % (-p))
        # We have polynomial convergence: use Richardson extrapolation
        # Need to use at least quad precision because a lot of cancellation
        # might occur in the extrapolation process
        prec2 = 4*prec
        one = MP_BASE(1) << prec2
        term = expr.subs(n, 0)
        term = (MP_BASE(term.p) << prec2) // term.q

        def summand(k, _term=[term]):
            if k:
                k = int(k)
                _term[0] *= MP_BASE(func1(k-1))
                _term[0] //= MP_BASE(func2(k-1))
            return make_mpf(from_man_exp(_term[0], -prec2))

        orig = mp.prec
        try:
            mp.prec = prec
            v = nsum(summand, [0, mpmath_inf], method='richardson')
        finally:
            mp.prec = orig
        return v._mpf_

def evalf_sum(expr, prec, options):
    func = expr.function
    limits = expr.limits
    if len(limits) != 1 or not isinstance(limits[0], tuple) or \
        len(limits[0]) != 3:
        raise NotImplementedError
    prec2 = prec+10
    try:
        n, a, b = limits[0]
        if b != S.Infinity or a != int(a):
            raise NotImplementedError
        # Use fast hypergeometric summation if possible
        v = hypsum(func, n, int(a), prec2)
        delta = prec - fastlog(v)
        if fastlog(v) < -10:
            v = hypsum(func, n, int(a), delta)
        return v, None, min(prec, delta), None
    except NotImplementedError:
        # Euler-Maclaurin summation for general series
        eps = C.Real(2.0)**(-prec)
        for i in range(1, 5):
            m = n = 2**i * prec
            s, err = expr.euler_maclaurin(m=m, n=n, eps=eps, \
                eval_integral=False)
            err = err.evalf()
            if err <= eps:
                break
        err = fastlog(evalf(abs(err), 20, options)[0])
        re, im, re_acc, im_acc = evalf(s, prec2, options)
        re_acc = max(re_acc, -err)
        im_acc = max(im_acc, -err)
        return re, im, re_acc, im_acc


#----------------------------------------------------------------------------#
#                                                                            #
#                            Symbolic interface                              #
#                                                                            #
#----------------------------------------------------------------------------#

def evalf_symbol(x, prec, options):
    val = options['subs'][x]
    if isinstance(val, mpf):
        if not val:
            return None, None, None, None
        return val._mpf_, None, prec, None
    else:
        if not '_cache' in options:
            options['_cache'] = {}
        cache = options['_cache']
        cached, cached_prec = cache.get(x.name, (None, MINUS_INF))
        if cached_prec >= prec:
            return cached
        v = evalf(sympify(val), prec, options)
        cache[x.name] = (v, prec)
        return v

evalf_table = None

def _create_evalf_table():
    global evalf_table
    evalf_table = {
    C.Symbol : evalf_symbol,
    C.Dummy : evalf_symbol,
    C.Real : lambda x, prec, options: (x._mpf_, None, prec, None),
    C.Rational : lambda x, prec, options: (from_rational(x.p, x.q, prec), None, prec, None),
    C.Integer : lambda x, prec, options: (from_int(x.p, prec), None, prec, None),
    C.Zero : lambda x, prec, options: (None, None, prec, None),
    C.One : lambda x, prec, options: (fone, None, prec, None),
    C.Half : lambda x, prec, options: (fhalf, None, prec, None),
    C.Pi : lambda x, prec, options: (mpf_pi(prec), None, prec, None),
    C.Exp1 : lambda x, prec, options: (mpf_e(prec), None, prec, None),
    C.ImaginaryUnit : lambda x, prec, options: (None, fone, None, prec),
    C.NegativeOne : lambda x, prec, options: (fnone, None, prec, None),

    C.exp : lambda x, prec, options: evalf_pow(C.Pow(S.Exp1, x.args[0],
        evaluate=False), prec, options),

    C.cos : evalf_trig,
    C.sin : evalf_trig,

    C.Add : evalf_add,
    C.Mul : evalf_mul,
    C.Pow : evalf_pow,

    C.log : evalf_log,
    C.atan : evalf_atan,
    C.abs : evalf_abs,

    C.re : evalf_re,
    C.im : evalf_im,
    C.floor : evalf_floor,
    C.ceiling : evalf_ceiling,

    C.Integral : evalf_integral,
    C.Sum : evalf_sum,
    C.Piecewise : evalf_piecewise,

    C.bernoulli : evalf_bernoulli,
    }

def evalf(x, prec, options):
    try:
        rf = evalf_table[x.func]
        r = rf(x, prec, options)
    except KeyError:
        #r = finalize_complex(x._eval_evalf(prec)._mpf_, fzero, prec)
        try:
            # Fall back to ordinary evalf if possible
            if 'subs' in options:
                x = x.subs(options['subs'])
            r = x._eval_evalf(prec)._mpf_, None, prec, None
        except AttributeError:
            raise NotImplementedError
    if options.get("verbose"):
        print "### input", x
        print "### output", to_str(r[0] or fzero, 50)
        print "### raw", r#r[0], r[2]
        print
    if options.get("chop"):
        r = chop_parts(r, prec)
    if options.get("strict"):
        check_target(x, r, prec)
    return r

def Basic_evalf(x, n=15, **options):
    """
    Evaluate the given formula to an accuracy of n digits.
    Optional keyword arguments:

        subs=<dict>
            Substitute numerical values for symbols, e.g.
            subs={x:3, y:1+pi}.

        maxprec=N
            Allow a maximum temporary working precision of N digits
            (default=100)

        chop=<bool>
            Replace tiny real or imaginary parts in subresults
            by exact zeros (default=False)

        strict=<bool>
            Raise PrecisionExhausted if any subresult fails to evaluate
            to full accuracy, given the available maxprec
            (default=False)

        quad=<str>
            Choose algorithm for numerical quadrature. By default,
            tanh-sinh quadrature is used. For oscillatory
            integrals on an infinite interval, try quad='osc'.

        verbose=<bool>
            Print debug information (default=False)

    """
    if not evalf_table:
        _create_evalf_table()
    prec = dps_to_prec(n)
    if 'maxprec' in options:
        options['maxprec'] = int(options['maxprec']*LG10)
    else:
        options['maxprec'] = max(prec, DEFAULT_MAXPREC)
    try:
        result = evalf(x, prec+4, options)
    except NotImplementedError:
        # Fall back to the ordinary evalf
        v = x._eval_evalf(prec)
        if v is None:
            return x
        try:
            # If the result is numerical, normalize it
            result = evalf(v, prec, options)
        except:
            # Probably contains symbols or unknown functions
            return v
    re, im, re_acc, im_acc = result
    if re:
        p = max(min(prec, re_acc), 1)
        #re = mpf_pos(re, p, round_nearest)
        re = C.Real._new(re, p)
    else:
        re = S.Zero
    if im:
        p = max(min(prec, im_acc), 1)
        #im = mpf_pos(im, p, round_nearest)
        im = C.Real._new(im, p)
        return re + im*S.ImaginaryUnit
    else:
        return re

Basic.evalf = Basic.n = Basic_evalf

def N(x, n=15, **options):
    """
    Calls x.evalf(n, **options).

    Both .evalf() and N() are equivalent, use the one that you like better.

    Example:
    >>> from sympy import Sum, Symbol, oo, N
    >>> from sympy.abc import k
    >>> Sum(1/k**k, (k, 1, oo))
    Sum(k**(-k), (k, 1, oo))
    >>> N(Sum(1/k**k, (k, 1, oo)), 4)
    1.291

    """
    return sympify(x).evalf(n, **options)
