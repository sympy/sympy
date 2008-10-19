"""
This module implements computation of hypergeometric and related
functions. In particular, it provides code for generic summation
of hypergeometric series. Optimized versions for various special
cases are also provided.
"""

import operator
import math

from settings import (\
    MP_ZERO, MP_ONE, round_fast, round_nearest
)

from libmpf import (\
    negative_rnd, bitcount, to_fixed, from_man_exp, to_int,
    fzero, fone, fnone, ftwo, finf, fninf, fnan,
    mpf_perturb, mpf_neg, mpf_shift, mpf_sub, mpf_mul, mpf_div,
    sqrt_fixed, mpf_sqrt, mpf_rdiv_int
)

from libelefun import (\
    mpf_pi, mpf_exp, pi_fixed
)

from libmpc import (\
    mpc_one, mpc_sub, mpc_mul_mpf, mpc_mul, mpc_neg, complex_int_pow
)

from gammazeta import int_fac

#-----------------------------------------------------------------------#
#                                                                       #
#                     Generic hypergeometric series                     #
#                                                                       #
#-----------------------------------------------------------------------#

"""
TODO:
  * By counting the number of multiplications vs divisions,
    the bit size of p can be kept around wp instead of growing
    it to n*wp for some (possibly large) n

  * Due to roundoff error, the series may fail to converge
    when x is negative and the convergence is slow.

"""

def hypsum_internal(ar, af, ac, br, bf, bc, xre, xim, prec, rnd):
    """
    Generic hypergeometric summation. This function computes:

            1   a_1 a_2 ...     1  (a_1 + 1) (a_2 + 1) ...  2
        1 + --  ----------- x + -- ----------------------- x  + ...
            1!  b_1 b_2 ...     2! (b_1 + 1) (b_2 + 1) ...

    The a_i and b_i sequences are separated by type:

    ar - list of a_i rationals [p,q]
    af - list of a_i mpf value tuples
    ac - list of a_i mpc value tuples
    br - list of b_i rationals [p,q]
    bf - list of b_i mpf value tuples
    bc - list of b_i mpc value tuples

    xre - real part of x, mpf value
    xim - imaginary part of x, mpf value (or None)

    Returns an mpc value is any complex parameter is given (or xim is not
    None); otherwise returns an mpf value.
    """

    have_float = af or bf
    have_complex = ac or bc

    # We want to mutate these in-place
    ar = [list(a) for a in ar]
    br = [list(b) for b in br]

    wp = prec + 25

    if xim is None:
        x = to_fixed(xre, wp)
        y = MP_ZERO
    else:
        have_complex = 1
        x = to_fixed(xre, wp)
        y = to_fixed(xim, wp)

    sre = pre = one = MP_ONE << wp
    sim = pim = MP_ZERO

    n = 1

    # Need to shift down by wp once for each fixed-point multiply
    # At minimum, we multiply by once by x each step
    shift = 1

    # Fixed-point real coefficients
    if have_float:
        len_af = len(af)
        len_bf = len(bf)
        range_af = range(len_af)
        range_bf = range(len_bf)
        for i in range_af: af[i] = to_fixed(af[i], wp)
        for i in range_bf: bf[i] = to_fixed(bf[i], wp)
        shift += len_af

    # Fixed-point complex coefficients
    if have_complex:
        len_ac = len(ac)
        len_bc = len(bc)
        range_ac = range(len_ac)
        range_bc = range(len_bc)
        for i in range_ac: ac[i] = [to_fixed(ac[i][0], wp), to_fixed(ac[i][1], wp)]
        for i in range_bc: bc[i] = [to_fixed(bc[i][0], wp), to_fixed(bc[i][1], wp)]
        shift += len_ac

    aqs = [a[1] for a in ar]
    bqs = [b[1] for b in br]
    aqprod = reduce(operator.mul, aqs, 1)
    bqprod = reduce(operator.mul, bqs, 1)

    assert shift >= 0

    while 1:
        # Integer and rational part of product
        mul = bqprod
        div = n * aqprod
        for ap, aq in ar: mul *= ap
        for bp, bq in br: div *= bp

        if have_complex:
            # Multiply by rational factors
            pre *= mul
            pim *= mul
            # Multiply by z
            pre, pim = pre*x - pim*y, pim*x + pre*y
            # Multiply by real factors
            for a in af:
                pre *= a
                pim *= a
            # Multiply by complex factors
            for are, aim in ac:
                pre, pim = pre*are - pim*aim, pim*are + pre*aim
            # Divide by rational factors
            pre //= div
            pim //= div
            # Divide by real factors
            for b in bf:
                pre = (pre << wp) // b
                pim = (pim << wp) // b
            # Divide by complex factors
            for bre, bim in bc:
                mag = bre*bre + bim*bim
                re = pre*bre + pim*bim
                im = pim*bre - pre*bim
                pre = (re << wp) // mag
                pim = (im << wp) // mag
        elif have_float:
            # Multiply and divide by real and rational factors, and x
            for a in af: pre *= a
            for b in bf:
                pre = (pre << wp) // b
            pre = (pre * (mul * x)) // div

        else:
            # Multiply and divide by rational factors and x
            pre = (pre * (mul * x)) // div

        pre >>= (wp*shift)
        sre += pre

        if have_complex:
            pim >>= (wp*shift)
            sim += pim
            if (-100 < pre < 100) and (-100 < pim < 100):
                break
        else:
            if -100 < pre < 100:
                break

        # Add 1 to all as and bs
        n += 1
        for ap_aq in ar: ap_aq[0] += ap_aq[1]
        for bp_bq in br: bp_bq[0] += bp_bq[1]
        if have_float:
            for i in range_af: af[i] += one
            for i in range_bf: bf[i] += one
        if have_complex:
            for i in range_ac: ac[i][0] += one
            for i in range_bc: bc[i][0] += one

    re = from_man_exp(sre, -wp, prec, rnd)
    if have_complex:
        im = from_man_exp(sim, -wp, prec, rnd)
        return (re, im)
    else:
        return re


#---------------------------------------------------------------------------#
#   Special-case implementation for rational parameters. These are          #
#   about 2x faster at low precision                                        #
#---------------------------------------------------------------------------#

def mpf_hyp0f1_rat((bp, bq), x, prec, rnd):
    wp = prec + 25
    x = to_fixed(x, wp)
    s = p = MP_ONE << wp
    n = 1
    while 1:
        p = (p * (bq*x) // (n*bp)) >> wp
        if -100 < p < 100:
            break
        s += p
        n += 1
        bp += bq
    return from_man_exp(s, -wp, prec, rnd)

def mpc_hyp0f1_rat((bp, bq), z, prec, rnd):
    wp = prec + 25
    zre, zim = z
    zre = to_fixed(zre, wp)
    zim = to_fixed(zim, wp)
    sre = pre = MP_ONE << wp
    sim = pim = MP_ZERO
    n = 1
    while 1:
        r1 = bq
        r2 = n*bp
        pre, pim = pre*zre - pim*zim, pim*zre + pre*zim
        pre = ((pre * r1) // r2) >> wp
        pim = ((pim * r1) // r2) >> wp
        if -100 < pre < 100 and -100 < pim < 100:
            break
        sre += pre
        sim += pim
        n += 1
        bp += bq
    re = from_man_exp(sre, -wp, prec, rnd)
    im = from_man_exp(sim, -wp, prec, rnd)
    return re, im

def mpf_hyp1f1_rat((ap, aq), (bp, bq), x, prec, rnd):
    wp = prec + 25
    x = to_fixed(x, wp)
    s = p = MP_ONE << wp
    n = 1
    while 1:
        p = (p * (ap*bq*x) // (n*aq*bp)) >> wp
        if -100 < p < 100:
            break
        s += p
        n += 1
        ap += aq
        bp += bq
    return from_man_exp(s, -wp, prec, rnd)

def mpc_hyp1f1_rat((ap, aq), (bp, bq), z, prec, rnd):
    wp = prec + 25
    zre, zim = z
    zre = to_fixed(zre, wp)
    zim = to_fixed(zim, wp)
    sre = pre = MP_ONE << wp
    sim = pim = MP_ZERO
    n = 1
    while 1:
        r1 = ap*bq
        r2 = n*aq*bp
        pre, pim = pre*zre - pim*zim, pim*zre + pre*zim
        pre = ((pre * r1) // r2) >> wp
        pim = ((pim * r1) // r2) >> wp
        if -100 < pre < 100 and -100 < pim < 100:
            break
        sre += pre
        sim += pim
        n += 1
        ap += aq
        bp += bq
    re = from_man_exp(sre, -wp, prec, rnd)
    im = from_man_exp(sim, -wp, prec, rnd)
    return re, im

def mpf_hyp2f1_rat((ap, aq), (bp, bq), (cp, cq), x, prec, rnd):
    wp = prec + 25
    x = to_fixed(x, wp)
    s = p = MP_ONE << wp
    n = 1
    while 1:
        p = (p * (ap*bp*cq*x) // (n*aq*bq*cp)) >> wp
        if -100 < p < 100:
            break
        s += p
        n += 1
        ap += aq
        bp += bq
        cp += cq
    return from_man_exp(s, -wp, prec, rnd)

def mpc_hyp2f1_rat((ap, aq), (bp, bq), (cp, cq), z, prec, rnd):
    wp = prec + 25
    zre, zim = z
    zre = to_fixed(zre, wp)
    zim = to_fixed(zim, wp)
    sre = pre = MP_ONE << wp
    sim = pim = MP_ZERO
    n = 1
    while 1:
        r1 = ap*bp*cq
        r2 = n*aq*bq*cp
        pre, pim = pre*zre - pim*zim, pim*zre + pre*zim
        pre = ((pre * r1) // r2) >> wp
        pim = ((pim * r1) // r2) >> wp
        if -100 < pre < 100 and -100 < pim < 100:
            break
        sre += pre
        sim += pim
        n += 1
        ap += aq
        bp += bq
        cp += cq
    re = from_man_exp(sre, -wp, prec, rnd)
    im = from_man_exp(sim, -wp, prec, rnd)
    return re, im



#-----------------------------------------------------------------------#
#                                                                       #
#                              Error functions                          #
#                                                                       #
#-----------------------------------------------------------------------#

# TODO: mpf_erf should call mpf_erfc when appropriate (currently
#    only the converse delegation is implemented)

def mpf_erf(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    if not man:
        if x == fzero: return fzero
        if x == finf: return fone
        if x== fninf: return fnone
        return fnan
    size = exp + bc
    lg = math.log
    # The approximation erf(x) = 1 is accurate to > x^2 * log(e,2) bits
    if size > 3 and 2*(size-1) + 0.528766 > lg(prec,2):
        if sign:
            return mpf_perturb(fnone, 0, prec, rnd)
        else:
            return mpf_perturb(fone, 1, prec, rnd)
    # erf(x) ~ 2*x/sqrt(pi) close to 0
    if size < -prec:
        # 2*x
        x = mpf_shift(x,1)
        c = mpf_sqrt(mpf_pi(prec+20), prec+20)
        # TODO: interval rounding
        return mpf_div(x, c, prec, rnd)
    wp = prec + abs(size) + 20
    # Taylor series for erf, fixed-point summation
    t = abs(to_fixed(x, wp))
    t2 = (t*t) >> wp
    s, term, k = t, 12345, 1
    while term:
        t = ((t * t2) >> wp) // k
        term = t // (2*k+1)
        if k & 1:
            s -= term
        else:
            s += term
        k += 1
    s = (s << (wp+1)) // sqrt_fixed(pi_fixed(wp), wp)
    if sign:
        s = -s
    return from_man_exp(s, -wp, wp, rnd)

def mpc_erf(z, prec, rnd=round_fast):
    re, im = z
    if im == fzero:
        return (mpf_erf(re, prec, rnd), fzero)
    wp = prec + 20
    z2 = mpc_mul(z, z, prec+20)
    v = mpc_hyp1f1_rat((1,2), (3,2), mpc_neg(z2), wp, rnd)
    sqrtpi = mpf_sqrt(mpf_pi(wp), wp)
    c = mpf_rdiv_int(2, sqrtpi, wp)
    c = mpc_mul_mpf(z, c, wp)
    return mpc_mul(c, v, prec, rnd)

# If possible, we use the asymptotic series for erfc.
# This is an alternating divergent asymptotic series, so
# the error is at most equal to the first omitted term.
# Here we check if the smallest term is small enough
# for a given x and precision
def erfc_check_series(x, prec):
    n = to_int(x)
    if n**2 * 1.44 > prec:
        return True
    return False

def mpf_erfc(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    if not man:
        if x == fzero: return fone
        if x == finf: return fzero
        if x == fninf: return ftwo
        return fnan
    wp = prec + 20
    mag = bc+exp
    # Preserve full accuracy when exponent grows huge
    wp += max(0, 2*mag)
    regular_erf = sign or mag < 2
    if regular_erf or not erfc_check_series(x, wp):
        if regular_erf:
            return mpf_sub(fone, mpf_erf(x, prec+10, negative_rnd[rnd]), prec, rnd)
        # 1-erf(x) ~ exp(-x^2), increase prec to deal with cancellation
        n = to_int(x)
        return mpf_sub(fone, mpf_erf(x, prec + int(n**2*1.44) + 10), prec, rnd)
    s = term = MP_ONE << wp
    term_prev = 0
    t = (2 * to_fixed(x, wp) ** 2) >> wp
    k = 1
    while 1:
        term = ((term * (2*k - 1)) << wp) // t
        if k > 4 and term > term_prev or not term:
            break
        if k & 1:
            s -= term
        else:
            s += term
        term_prev = term
        #print k, to_str(from_man_exp(term, -wp, 50), 10)
        k += 1
    s = (s << wp) // sqrt_fixed(pi_fixed(wp), wp)
    s = from_man_exp(s, -wp, wp)
    z = mpf_exp(mpf_neg(mpf_mul(x,x,wp),wp),wp)
    y = mpf_div(mpf_mul(z, s, wp), x, prec, rnd)
    return y

def mpc_erfc(z, prec, rnd=round_fast):
    real, imag = z
    if not imag:
        return (mpf_erfc(real, prec, rnd), fzero)
    # XXX: cancellation
    return mpc_sub(mpc_one, mpc_erf(z, prec+20, rnd), prec, rnd)

#-----------------------------------------------------------------------#
#                                                                       #
#                             Bessel functions                          #
#                                                                       #
#-----------------------------------------------------------------------#

# A Bessel function of the first kind of integer order, J_n(x), is
# given by the power series

#             oo
#             ___         k         2 k + n
#            \        (-1)     / x \
#    J_n(x) = )    ----------- | - |
#            /___  k! (k + n)! \ 2 /
#            k = 0

# Simplifying the quotient between two successive terms gives the
# ratio x^2 / (-4*k*(k+n)). Hence, we only need one full-precision
# multiplication and one division by a small integer per term.
# The complex version is very similar, the only difference being
# that the multiplication is actually 4 multiplies.

# In the general case, we have
# J_v(x) = (x/2)**v / v! * 0F1(v+1, (-1/4)*z**2)

# TODO: for extremely large x, we could use an asymptotic
# trigonometric approximation.

# TODO: recompute at higher precision if the fixed-point mantissa
# is very small

def mpf_besseljn(n, x, prec):
    negate = n < 0 and n & 1
    n = abs(n)
    origprec = prec
    prec += 20 + bitcount(abs(n))
    x = to_fixed(x, prec)
    x2 = (x**2) >> prec
    if not n:
        s = t = MP_ONE << prec
    else:
        s = t = (x**n // int_fac(n)) >> ((n-1)*prec + n)
    k = 1
    while t:
        t = ((t * x2) // (-4*k*(k+n))) >> prec
        s += t
        k += 1
    if negate:
        s = -s
    return from_man_exp(s, -prec, origprec, round_nearest)

def mpc_besseljn(n, z, prec):
    negate = n < 0 and n & 1
    n = abs(n)
    origprec = prec
    prec += 20 + bitcount(abs(n))
    zre, zim = z
    zre = to_fixed(zre, prec)
    zim = to_fixed(zim, prec)
    z2re = (zre**2 - zim**2) >> prec
    z2im = (zre*zim) >> (prec-1)
    if not n:
        sre = tre = MP_ONE << prec
        sim = tim = MP_ZERO
    else:
        re, im = complex_int_pow(zre, zim, n)
        sre = tre = (re // int_fac(n)) >> ((n-1)*prec + n)
        sim = tim = (im // int_fac(n)) >> ((n-1)*prec + n)
    k = 1
    while abs(tre) + abs(tim) > 3:
        p = -4*k*(k+n)
        tre, tim = tre*z2re - tim*z2im, tim*z2re + tre*z2im
        tre = (tre // p) >> prec
        tim = (tim // p) >> prec
        sre += tre
        sim += tim
        k += 1
    if negate:
        sre = -sre
        sim = -sim
    re = from_man_exp(sre, -prec, origprec, round_nearest)
    im = from_man_exp(sim, -prec, origprec, round_nearest)
    return (re, im)
