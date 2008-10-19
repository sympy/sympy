"""
Low-level functions for complex arithmetic.
"""

from settings import (\
    MP_BASE, MP_ZERO, MP_ONE, MP_TWO,
    round_floor, round_ceiling, round_down, round_up,
    round_nearest, round_fast
)

from libmpf import (\
    bctable, normalize, reciprocal_rnd, rshift, lshift, giant_steps,
    to_str, to_fixed, from_man_exp, from_float, from_int, to_int,
    fzero, fone, ftwo, fhalf, finf, fninf, fnan,
    mpf_abs, mpf_pos, mpf_neg, mpf_add, mpf_sub, mpf_mul,
    mpf_div, mpf_mul_int, mpf_shift, mpf_sqrt, mpf_hypot,
    mpf_rdiv_int, mpf_floor, mpf_ceil
)

from libelefun import (\
    mpf_pi, mpf_exp, mpf_log, cos_sin, cosh_sinh, mpf_tan,
    mpf_atan, mpf_atan2, mpf_cosh, mpf_sinh, mpf_tanh,
    mpf_asin, mpf_acos, mpf_acosh
)


# An mpc value is a (real, imag) tuple
mpc_one = fone, fzero
mpc_zero = fzero, fzero
mpc_two = ftwo, fzero
mpc_half = (fhalf, fzero)

def complex_to_str(re, im, dps):
    rs = to_str(re, dps)
    if im[0]:
        return rs + " - " + to_str(mpf_neg(im), dps) + "j"
    else:
        return rs + " + " + to_str(im, dps) + "j"

def mpc_add((a, b), (c, d), prec, rnd=round_fast):
    return mpf_add(a, c, prec, rnd), mpf_add(b, d, prec, rnd)

def mpc_add_mpf((a, b), p, prec, rnd=round_fast):
      return mpf_add(a, p, prec, rnd), b

def mpc_sub((a, b), (c, d), prec, rnd=round_fast):
    return mpf_sub(a, c, prec, rnd), mpf_sub(b, d, prec, rnd)

def mpc_sub_mpf((a, b), p, prec, rnd=round_fast):
    return mpf_sub(a, p, prec, rnd), b

def mpc_pos((a, b), prec, rnd=round_fast):
    return mpf_pos(a, prec, rnd), mpf_pos(b, prec, rnd)

def mpc_neg((a, b), prec=None, rnd=round_fast):
    return mpf_neg(a, prec, rnd), mpf_neg(b, prec, rnd)

def mpc_shift((a, b), n):
    return mpf_shift(a, n), mpf_shift(b, n)

def mpc_abs((a, b), prec, rnd=round_fast):
    """Absolute value of a complex number, |a+bi|.
    Returns an mpf value."""
    return mpf_hypot(a, b, prec, rnd)

def mpc_arg((a, b), prec, rnd=round_fast):
    """Argument of a complex number. Returns an mpf value."""
    return mpf_atan2(b, a, prec, rnd)

def mpc_floor((a, b), prec, rnd=round_fast):
    return mpf_floor(a, prec, rnd), mpf_floor(b, prec, rnd)

def mpc_ceil((a, b), prec, rnd=round_fast):
    return mpf_ceil(a, prec, rnd), mpf_ceil(b, prec, rnd)

def mpc_mul((a, b), (c, d), prec, rnd=round_fast):
    """Complex multiplication.

    Returns the real and imaginary part of (a+bi)*(c+di), rounded to
    the specified precision. The rounding mode applies to the real and
    imaginary parts separately."""

    asign, aman, aexp, abc = a
    bsign, bman, bexp, bbc = b
    csign, cman, cexp, cbc = c
    dsign, dman, dexp, dbc = d

    if 0 in (aman, bman, cman, dman):
        # TODO: use a better strategy for complex infs
        if ((not aman) and aexp) or ((not bman) and bexp) or \
            ((not cman) and cexp) or ((not dman) and dexp):
            wp = prec + 10
            re = mpf_sub(mpf_mul(a,c), mpf_mul(b,d), prec, rnd)
            im = mpf_add(mpf_mul(a,d), mpf_mul(b,c), prec, rnd)
            return re, im
        # bi * c = bci
        # bi * di = -bd
        # bi * (c + di) = -bd + bci
        if not aman:
            if not dman: return fzero, mpf_mul(b, c, prec, rnd)
            if not cman: return mpf_mul(mpf_neg(b), d, prec, rnd), fzero
            return mpf_mul(mpf_neg(b), d, prec, rnd), mpf_mul(b, c, prec, rnd)
        # a * c = ac
        # a * di = adi
        # a * (c + di) = ac + adi
        if not bman:
            if not dman: return mpf_mul(a, c, prec, rnd), fzero
            if not cman: return fzero, mpf_mul(a, d, prec, rnd)
            return mpf_mul(a, c, prec, rnd), mpf_mul(a, d, prec, rnd)
        # (a + bi) * c
        if not dman:
            return mpf_mul(a, c, prec, rnd), mpf_mul(b, c, prec, rnd)
        # (a + bi) * di = -bd + adi
        return mpf_mul(mpf_neg(b), d, prec, rnd), mpf_mul(a, d, prec, rnd)

    # Avoid normalizing the temporary products
    bct = bctable

    psign = asign ^ csign
    pman = aman * cman
    pexp = aexp + cexp
    pbc = abc + cbc - 4
    if pbc < 4: pbc = bct[int(pman)]
    else:       pbc += bct[int(pman>>pbc)]
    p = psign, pman, pexp, pbc

    qsign = (bsign ^ dsign) ^ 1
    qman = bman * dman
    qexp = bexp + dexp
    qbc = bbc + dbc - 4
    if qbc < 4: qbc = bct[int(qman)]
    else:       qbc += bct[int(qman>>qbc)]
    q = qsign, qman, qexp, qbc

    rsign = bsign ^ csign
    rman = bman * cman
    rexp = bexp + cexp
    rbc = bbc + cbc - 4
    if rbc < 4: rbc = bct[int(rman)]
    else:       rbc += bct[int(rman>>rbc)]
    r = rsign, rman, rexp, rbc

    ssign = asign ^ dsign
    sman = aman * dman
    sexp = aexp + dexp
    sbc = abc + dbc - 4
    if sbc < 4: sbc = bct[int(sman)]
    else:       sbc += bct[int(sman>>sbc)]
    s = ssign, sman, sexp, sbc

    return mpf_add(p, q, prec, rnd), mpf_add(r, s, prec, rnd)


def mpc_mul_mpf((a, b), p, prec, rnd=round_fast):
    re = mpf_mul(a, p, prec, rnd)
    im = mpf_mul(b, p, prec, rnd)
    return re, im

def mpc_mul_int((a, b), n, prec, rnd=round_fast):
    re = mpf_mul_int(a, n, prec, rnd)
    im = mpf_mul_int(b, n, prec, rnd)
    return re, im

def mpc_div((a, b), (c, d), prec, rnd=round_fast):
    wp = prec + 10
    # mag = c*c + d*d
    mag = mpf_add(mpf_mul(c, c), mpf_mul(d, d), wp)
    # (a*c+b*d)/mag, (b*c-a*d)/mag
    t = mpf_add(mpf_mul(a,c), mpf_mul(b,d), wp)
    u = mpf_sub(mpf_mul(b,c), mpf_mul(a,d), wp)
    return mpf_div(t,mag,prec,rnd), mpf_div(u,mag,prec,rnd)

def mpc_div_mpf((a, b), p, prec, rnd=round_fast):
    re = mpf_div(a, p, prec, rnd)
    im = mpf_div(b, p, prec, rnd)
    return re, im

def complex_int_pow(a, b, n):
    """Complex integer power: computes (a+b*I)**n exactly for
    nonnegative n (a and b must be Python ints)."""
    wre = 1
    wim = 0
    while n:
        if n & 1:
            wre, wim = wre*a - wim*b, wim*a + wre*b
            n -= 1
        a, b = a*a - b*b, 2*a*b
        n //= 2
    return wre, wim

def mpc_pow(z, w, prec, rnd=round_fast):
    if w[1] == fzero:
        return mpc_pow_mpf(z, w[0], prec, rnd)
    return mpc_exp(mpc_mul(mpc_log(z, prec+10), w, prec+10), prec, rnd)

def mpc_pow_mpf(z, p, prec, rnd=round_fast):
    psign, pman, pexp, pbc = p
    if pexp >= 0:
        return mpc_pow_int(z, (-1)**psign * (pman<<pexp), prec, rnd)
    if pexp == -1:
        sqrtz = mpc_sqrt(z, prec+10)
        return mpc_pow_int(sqrtz, (-1)**psign * pman, prec, rnd)
    return mpc_exp(mpc_mul_mpf(mpc_log(z, prec+10), p, prec+10), prec, rnd)

def mpc_pow_int(z, n, prec, rnd=round_fast):
    if n == 0: return mpc_one
    if n == 1: return mpc_pos(z, prec, rnd)
    if n == 2: return mpc_mul(z, z, prec, rnd)
    if n == -1: return mpc_div(mpc_one, z, prec, rnd)
    if n < 0: return mpc_div(mpc_one, mpc_pow_int(z, -n, prec+4), prec, rnd)
    a, b = z
    asign, aman, aexp, abc = a
    bsign, bman, bexp, bbc = b
    if asign: aman = -aman
    if bsign: bman = -bman
    de = aexp - bexp
    abs_de = abs(de)
    exact_size = n*(abs_de + max(abc, bbc))
    if exact_size < 10000:
        if de > 0:
            aman <<= de
            aexp = bexp
        else:
            bman <<= (-de)
            bexp = aexp
        re, im = complex_int_pow(aman, bman, n)
        re = from_man_exp(re, int(n*aexp), prec, rnd)
        im = from_man_exp(im, int(n*bexp), prec, rnd)
        return re, im
    return mpc_exp(mpc_mul_int(mpc_log(z, prec+10), n, prec+10), prec, rnd)

def mpc_sqrt((a, b), prec, rnd=round_fast):
    """Complex square root (principal branch).

    We have sqrt(a+bi) = sqrt((r+a)/2) + b/sqrt(2*(r+a))*i where
    r = abs(a+bi), when a+bi is not a negative real number."""
    if a == b == fzero:
        return (a, b)
    # When a+bi is a negative real number, we get a real sqrt times i
    if b == fzero:
        if a[0]:
            im = mpf_sqrt(mpf_neg(a), prec, rnd)
            return (fzero, im)
        else:
            re = mpf_sqrt(a, prec, rnd)
            return (re, fzero)
    wp = prec+20
    if not a[0]:                               # case a positive
        t  = mpf_add(mpc_abs((a, b), wp), a, wp)  # t = abs(a+bi) + a
        u = mpf_shift(t, -1)                      # u = t/2
        re = mpf_sqrt(u, prec, rnd)               # re = sqrt(u)
        v = mpf_shift(t, 1)                       # v = 2*t
        w  = mpf_sqrt(v, wp)                      # w = sqrt(v)
        im = mpf_div(b, w, prec, rnd)             # im = b / w
    else:                                      # case a negative
        t = mpf_sub(mpc_abs((a, b), wp), a, wp)   # t = abs(a+bi) - a
        u = mpf_shift(t, -1)                      # u = t/2
        im = mpf_sqrt(u, prec, rnd)               # im = sqrt(u)
        v = mpf_shift(t, 1)                       # v = 2*t
        w  = mpf_sqrt(v, wp)                      # w = sqrt(v)
        re = mpf_div(b, w, prec, rnd)             # re = b/w
        if b[0]:
            re = mpf_neg(re)
            im = mpf_neg(im)
    return re, im

def mpc_nthroot_fixed(a, b, n, prec):
    # a, b signed integers at fixed precision prec
    start = 50
    a1 = int(rshift(a, prec - n*start))
    b1 = int(rshift(b, prec - n*start))
    try:
        r = (a1 + 1j * b1)**(1.0/n)
        re = r.real
        im = r.imag
        # XXX: workaround bug in gmpy
        if abs(re) < 0.1: re = 0
        if abs(im) < 0.1: im = 0
        re = MP_BASE(re)
        im = MP_BASE(im)
    except OverflowError:
        a1 = from_int(a1, start)
        b1 = from_int(b1, start)
        fn = from_int(n)
        nth = mpf_rdiv_int(1, fn, start)
        re, im = mpc_pow((a1, b1), (nth, fzero), start)
        re = to_int(re)
        im = to_int(im)
    extra = 10
    prevp = start
    extra1 = n
    for p in giant_steps(start, prec+extra):
        # this is slow for large n, unlike int_pow_fixed
        re2, im2 = complex_int_pow(re, im, n-1)
        re2 = rshift(re2, (n-1)*prevp - p - extra1)
        im2 = rshift(im2, (n-1)*prevp - p - extra1)
        r4 = (re2*re2 + im2*im2) >> (p + extra1)
        ap = rshift(a, prec - p)
        bp = rshift(b, prec - p)
        rec = (ap * re2 + bp * im2) >> p
        imc = (-ap * im2 + bp * re2) >> p
        reb = (rec << p) // r4
        imb = (imc << p) // r4
        re = (reb + (n-1)*lshift(re, p-prevp))//n
        im = (imb + (n-1)*lshift(im, p-prevp))//n
        prevp = p
    return re, im

def mpc_nthroot((a, b), n, prec, rnd=round_fast):
    """
    Complex n-th root.

    Use Newton method as in the real case when it is faster,
    otherwise use z**(1/n)
    """

    if a[0] == 0 and b == fzero:
        re = mpf_nthroot(a, n, prec, rnd)
        return (re, fzero)
    if n < 2:
        if n == 0:
            return mpc_one
        if n == 1:
            return mpc_pos((a, b), prec, rnd)
        if n == -1:
            return mpc_div(mpc_one, (a, b), prec, rnd)
        inverse = mpc_nthroot((a, b), -n, prec+5, reciprocal_rnd[rnd])
        return mpc_div(mpc_one, inverse, prec, rnd)
    if n > 20:
        fn = from_int(n)
        prec2 = prec+10
        nth = mpf_rdiv_int(1, fn, prec2)
        re, im = mpc_pow((a, b), (nth, fzero), prec2, rnd)
        re = normalize(re[0], re[1], re[2], re[3], prec, rnd)
        im = normalize(im[0], im[1], im[2], im[3], prec, rnd)
        return re, im
    prec2 = int(1.2 * (prec + 10))
    asign, aman, aexp, abc = a
    bsign, bman, bexp, bbc = b
    af = to_fixed(a, prec2)
    bf = to_fixed(b, prec2)
    re, im = mpc_nthroot_fixed(af, bf, n, prec2)
    extra = 10
    re = from_man_exp(re, -prec2-extra, prec2, rnd)
    im = from_man_exp(im, -prec2-extra, prec2, rnd)
    return re, im

def mpc_cbrt((a, b), prec, rnd=round_fast):
    """
    Complex cubic root.
    """
    return mpc_nthroot((a, b), 3, prec, rnd)

def mpc_exp((a, b), prec, rnd=round_fast):
    """
    Complex exponential function.

    We use the direct formula exp(a+bi) = exp(a) * (cos(b) + sin(b)*i)
    for the computation. This formula is very nice because it is
    pewrectly stable; since we just do real multiplications, the only
    numerical errors that can crewp in are single-ulp rnd errors.

    The formula is efficient since mpmath's real exp is quite fast and
    since we can compute cos and sin simultaneously.

    It is no problem if a and b are large; if the implementations of
    exp/cos/sin are accurate and efficient for all real numbers, then
    so is this function for all complex numbers.
    """
    if a == fzero:
        return cos_sin(b, prec, rnd)
    mag = mpf_exp(a, prec+4, rnd)
    c, s = cos_sin(b, prec+4, rnd)
    re = mpf_mul(mag, c, prec, rnd)
    im = mpf_mul(mag, s, prec, rnd)
    return re, im

def mpc_log(z, prec, rnd=round_fast):
    return mpf_log(mpc_abs(z, prec, rnd), prec, rnd), mpc_arg(z, prec, rnd)

def mpc_cos((a, b), prec, rnd=round_fast):
    """Complex cosine. The formula used is cos(a+bi) = cos(a)*cosh(b) -
    sin(a)*sinh(b)*i.

    The same comments apply as for the complex exp: only real
    multiplications are pewrormed, so no cancellation errors are
    possible. The formula is also efficient since we can compute both
    pairs (cos, sin) and (cosh, sinh) in single stwps."""
    if a == fzero:
        return mpf_cosh(b, prec, rnd), fzero
    wp = prec + 6
    c, s = cos_sin(a, wp)
    ch, sh = cosh_sinh(b, wp)
    re = mpf_mul(c, ch, prec, rnd)
    im = mpf_mul(s, sh, prec, rnd)
    return re, mpf_neg(im)

def mpc_sin((a, b), prec, rnd=round_fast):
    """Complex sine. We have sin(a+bi) = sin(a)*cosh(b) +
    cos(a)*sinh(b)*i. See the docstring for mpc_cos for additional
    comments."""
    if a == fzero:
        return fzero, mpf_sinh(b, prec, rnd)
    wp = prec + 6
    c, s = cos_sin(a, wp)
    ch, sh = cosh_sinh(b, wp)
    re = mpf_mul(s, ch, prec, rnd)
    im = mpf_mul(c, sh, prec, rnd)
    return re, im

def mpc_tan(z, prec, rnd=round_fast):
    """Complex tangent. Computed as tan(a+bi) = sin(2a)/M + sinh(2b)/M*i
    where M = cos(2a) + cosh(2b)."""
    a, b = z
    asign, aman, aexp, abc = a
    bsign, bman, bexp, bbc = b
    if b == fzero: return mpf_tan(a, prec, rnd), fzero
    if a == fzero: return fzero, mpf_tanh(b, prec, rnd)
    wp = prec + 15
    a = mpf_shift(a, 1)
    b = mpf_shift(b, 1)
    c, s = cos_sin(a, wp)
    ch, sh = cosh_sinh(b, wp)
    # TODO: handle cancellation when c ~=  -1 and ch ~= 1
    mag = mpf_add(c, ch, wp)
    re = mpf_div(s, mag, prec, rnd)
    im = mpf_div(sh, mag, prec, rnd)
    return re, im

def mpc_cosh((a, b), prec, rnd=round_fast):
    """Complex hyperbolic cosine. Computed as cosh(z) = cos(z*i)."""
    return mpc_cos((b, mpf_neg(a)), prec, rnd)

def mpc_sinh((a, b), prec, rnd=round_fast):
    """Complex hyperbolic sine. Computed as sinh(z) = -i*sin(z*i)."""
    b, a = mpc_sin((b, a), prec, rnd)
    return a, b

def mpc_tanh((a, b), prec, rnd=round_fast):
    """Complex hyperbolic tangent. Computed as tanh(z) = -i*tan(z*i)."""
    b, a = mpc_tan((b, a), prec, rnd)
    return a, b

# TODO: avoid loss of accuracy
def mpc_atan((a, b), prec, rnd=round_fast):
    # atan(z) = (I/2)*(log(1-I*z) - log(1+I*z))
    # x = 1-I*z = 1 + b - I*a
    # y = 1+I*z = 1 - b + I*a
    wp = prec + 15
    x = mpf_add(fone, b, wp), mpf_neg(a)
    y = mpf_sub(fone, b, wp), a
    l1 = mpc_log(x, wp)
    l2 = mpc_log(y, wp)
    a, b = mpc_sub(l1, l2, prec, rnd)
    # (I/2) * (a+b*I) = (-b/2 + a/2*I)
    return mpf_neg(mpf_shift(b,-1)), mpf_shift(a,-1)

beta_crossover = from_float(0.6417)
alpha_crossover = from_float(1.5)

def acos_asin(z, prec, rnd, n):
    """ complex acos for n = 0, asin for n = 1
    The algorithm is described in
    T.E. Hull, T.F. Fairgrieve and P.T.P. Tang
    'Implementing the Complex Arcsine and Arcosine Functions
    using Exception Handling',
    ACM Trans. on Math. Software Vol. 23 (1997), p299
    The complex acos and asin can be defined as
    acos(z) = acos(beta) - I*sign(a)* log(alpha + sqrt(alpha**2 -1))
    asin(z) = asin(beta) + I*sign(a)* log(alpha + sqrt(alpha**2 -1))
    where z = a + I*b
    alpha = (1/2)*(r + s); beta = (1/2)*(r - s) = a/alpha
    r = sqrt((a+1)**2 + y**2); s = sqrt((a-1)**2 + y**2)
    These expressions are rewritten in different ways in different
    regions, delimited by two crossovers alpha_crossover and beta_crossover,
    and by abs(a) <= 1, in order to improve the numerical accuracy.
    """
    a, b = z
    wp = prec + 10
    # special cases with real argument
    if b == fzero:
        am = mpf_sub(fone, mpf_abs(a), wp)
        # case abs(a) <= 1
        if not am[0]:
            if n == 0:
                return mpf_acos(a, prec, rnd), fzero
            else:
                return mpf_asin(a, prec, rnd), fzero
        # cases abs(a) > 1
        else:
            # case a < -1
            if a[0]:
                pi = mpf_pi(prec, rnd)
                c = mpf_acosh(mpf_neg(a), prec, rnd)
                if n == 0:
                    return pi, mpf_neg(c)
                else:
                    return mpf_neg(mpf_shift(pi, -1)), c
            # case a > 1
            else:
                c = mpf_acosh(a, prec, rnd)
                if n == 0:
                    return fzero, c
                else:
                    pi = mpf_pi(prec, rnd)
                    return mpf_shift(pi, -1), mpf_neg(c)
    asign = bsign = 0
    if a[0]:
        a = mpf_neg(a)
        asign = 1
    if b[0]:
        b = mpf_neg(b)
        bsign = 1
    am = mpf_sub(fone, a, wp)
    ap = mpf_add(fone, a, wp)
    r = mpf_hypot(ap, b, wp)
    s = mpf_hypot(am, b, wp)
    alpha = mpf_shift(mpf_add(r, s, wp), -1)
    beta = mpf_div(a, alpha, wp)
    b2 = mpf_mul(b,b, wp)
    # case beta <= beta_crossover
    if not mpf_sub(beta_crossover, beta, wp)[0]:
        if n == 0:
            re = mpf_acos(beta, wp)
        else:
            re = mpf_asin(beta, wp)
    else:
        # to compute the real part in this region use the identity
        # asin(beta) = atan(beta/sqrt(1-beta**2))
        # beta/sqrt(1-beta**2) = (alpha + a) * (alpha - a)
        # alpha + a is numerically accurate; alpha - a can have
        # cancellations leading to numerical inaccuracies, so rewrite
        # it in differente ways according to the region
        Ax = mpf_add(alpha, a, wp)
        # case a <= 1
        if not am[0]:
            # c = b*b/(r + (a+1)); d = (s + (1-a))
            # alpha - a = (1/2)*(c + d)
            # case n=0: re = atan(sqrt((1/2) * Ax * (c + d))/a)
            # case n=1: re = atan(a/sqrt((1/2) * Ax * (c + d)))
            c = mpf_div(b2, mpf_add(r, ap, wp), wp)
            d = mpf_add(s, am, wp)
            re = mpf_shift(mpf_mul(Ax, mpf_add(c, d, wp), wp), -1)
            if n == 0:
                re = mpf_atan(mpf_div(mpf_sqrt(re, wp), a, wp), wp)
            else:
                re = mpf_atan(mpf_div(a, mpf_sqrt(re, wp), wp), wp)
        else:
            # c = Ax/(r + (a+1)); d = Ax/(s - (1-a))
            # alpha - a = (1/2)*(c + d)
            # case n = 0: re = atan(b*sqrt(c + d)/2/a)
            # case n = 1: re = atan(a/(b*sqrt(c + d)/2)
            c = mpf_div(Ax, mpf_add(r, ap, wp), wp)
            d = mpf_div(Ax, mpf_sub(s, am, wp), wp)
            re = mpf_shift(mpf_add(c, d, wp), -1)
            re = mpf_mul(b, mpf_sqrt(re, wp), wp)
            if n == 0:
                re = mpf_atan(mpf_div(re, a, wp), wp)
            else:
                re = mpf_atan(mpf_div(a, re, wp), wp)
    # to compute alpha + sqrt(alpha**2 - 1), if alpha <= alpha_crossover
    # replace it with 1 + Am1 + sqrt(Am1*(alpha+1)))
    # where Am1 = alpha -1
    # if alpha <= alpha_crossover:
    if not mpf_sub(alpha_crossover, alpha, wp)[0]:
        c1 = mpf_div(b2, mpf_add(r, ap, wp), wp)
        # case a < 1
        if mpf_neg(am)[0]:
            # Am1 = (1/2) * (b*b/(r + (a+1)) + b*b/(s + (1-a))
            c2 = mpf_add(s, am, wp)
            c2 = mpf_div(b2, c2, wp)
            Am1 = mpf_shift(mpf_add(c1, c2, wp), -1)
        else:
            # Am1 = (1/2) * (b*b/(r + (a+1)) + (s - (1-a)))
            c2 = mpf_sub(s, am, wp)
            Am1 = mpf_shift(mpf_add(c1, c2, wp), -1)
        # im = log(1 + Am1 + sqrt(Am1*(alpha+1)))
        im = mpf_mul(Am1, mpf_add(alpha, fone, wp), wp)
        im = mpf_log(mpf_add(fone, mpf_add(Am1, mpf_sqrt(im, wp), wp), wp), wp)
    else:
        # im = log(alpha + sqrt(alpha*alpha - 1))
        im = mpf_sqrt(mpf_sub(mpf_mul(alpha, alpha, wp), fone, wp), wp)
        im = mpf_log(mpf_add(alpha, im, wp), wp)
    if asign:
        if n == 0:
            re = mpf_sub(mpf_pi(wp), re, wp)
        else:
            re = mpf_neg(re)
    if not bsign and n == 0:
        im = mpf_neg(im)
    if bsign and n == 1:
        im = mpf_neg(im)
    re = normalize(re[0], re[1], re[2], re[3], prec, rnd)
    im = normalize(im[0], im[1], im[2], im[3], prec, rnd)
    return re, im

def mpc_acos(z, prec, rnd=round_fast):
    return acos_asin(z, prec, rnd, 0)

def mpc_asin(z, prec, rnd=round_fast):
    return acos_asin(z, prec, rnd, 1)

def mpc_asinh(z, prec, rnd=round_fast):
    # asinh(z) = I * asin(-I z)
    a, b = z
    a, b =  mpc_asin((b, mpf_neg(a)), prec, rnd)
    return mpf_neg(b), a

def mpc_acosh(z, prec, rnd=round_fast):
    # acosh(z) = -I * acos(z)   for Im(acos(z)) <= 0
    #            +I * acos(z)   otherwise
    a, b = mpc_acos(z, prec, rnd)
    if b[0] or b == fzero:
        return mpf_neg(b), a
    else:
        return b, mpf_neg(a)

def mpc_atanh(z, prec, rnd=round_fast):
    # atanh(z) = (log(1+z)-log(1-z))/2
    wp = prec + 15
    a = mpc_add(z, mpc_one, wp)
    b = mpc_sub(mpc_one, z, wp)
    a = mpc_log(a, wp)
    b = mpc_log(b, wp)
    return mpc_shift(mpc_sub(a, b, wp), -1)
