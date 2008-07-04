"""
Low-level functions for complex arithmetic.
"""

from lib import *

# An mpc value is a (real, imag) tuple
mpc_one = fone, fzero
mpc_zero = fzero, fzero
mpc_two = ftwo, fzero

def complex_to_str(re, im, dps):
    rs = to_str(re, dps)
    if im[0]:
        return rs + " - " + to_str(fneg(im), dps) + "j"
    else:
        return rs + " + " + to_str(im, dps) + "j"

def mpc_add((a, b), (c, d), prec, rnd=round_fast):
    return fadd(a, c, prec, rnd), fadd(b, d, prec, rnd)

def mpc_add_mpf((a, b), p, prec, rnd=round_fast):
      return fadd(a, p, prec, rnd), b

def mpc_sub((a, b), (c, d), prec, rnd=round_fast):
    return fsub(a, c, prec, rnd), fsub(b, d, prec, rnd)

def mpc_sub_mpf((a, b), p, prec, rnd=round_fast):
    return fsub(a, p, prec, rnd), b

def mpc_pos((a, b), prec, rnd=round_fast):
    return fpos(a, prec, rnd), fpos(b, prec, rnd)

def mpc_shift((a, b), n):
    return fshift(a, n), fshift(b, n)

def mpc_abs((a, b), prec, rnd=round_fast):
    """Absolute value of a complex number, |a+bi|.
    Returns an mpf value."""
    return fhypot(a, b, prec, rnd)

def mpc_arg((a, b), prec, rnd=round_fast):
    """Argument of a complex number. Returns an mpf value."""
    return fatan2(b, a, prec, rnd)

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
            re = fsub(fmul(a,c), fmul(b,d), prec, rnd)
            im = fadd(fmul(a,d), fmul(b,c), prec, rnd)
            return re, im
        # bi * c = bci
        # bi * di = -bd
        # bi * (c + di) = -bd + bci
        if not aman:
            if not dman: return fzero, fmul(b, c, prec, rnd)
            if not cman: return fmul(fneg(b), d, prec, rnd), fzero
            return fmul(fneg(b), d, prec, rnd), fmul(b, c, prec, rnd)
        # a * c = ac
        # a * di = adi
        # a * (c + di) = ac + adi
        if not bman:
            if not dman: return fmul(a, c, prec, rnd), fzero
            if not cman: return fzero, fmul(a, d, prec, rnd)
            return fmul(a, c, prec, rnd), fmul(a, d, prec, rnd)
        # (a + bi) * c
        if not dman:
            return fmul(a, c, prec, rnd), fmul(b, c, prec, rnd)
        # (a + bi) * di = -bd + adi
        return fmul(fneg(b), d, prec, rnd), fmul(a, d, prec, rnd)

    # Avoid normalizing the temporary products
    bct = bctable

    psign = asign ^ csign
    pman = aman * cman
    pexp = aexp + cexp
    pbc = abc + cbc - 4
    if pbc < 4: pbc = bct[pman]
    else:       pbc += bct[pman>>pbc]
    p = psign, pman, pexp, pbc

    qsign = (bsign ^ dsign) ^ 1
    qman = bman * dman
    qexp = bexp + dexp
    qbc = bbc + dbc - 4
    if qbc < 4: qbc = bct[qman]
    else:       qbc += bct[qman>>qbc]
    q = qsign, qman, qexp, qbc

    rsign = bsign ^ csign
    rman = bman * cman
    rexp = bexp + cexp
    rbc = bbc + cbc - 4
    if rbc < 4: rbc = bct[rman]
    else:       rbc += bct[rman>>rbc]
    r = rsign, rman, rexp, rbc

    ssign = asign ^ dsign
    sman = aman * dman
    sexp = aexp + dexp
    sbc = abc + dbc - 4
    if sbc < 4: sbc = bct[sman]
    else:       sbc += bct[sman>>sbc]
    s = ssign, sman, sexp, sbc

    return fadd(p, q, prec, rnd), fadd(r, s, prec, rnd)


def mpc_mul_mpf((a, b), p, prec, rnd=round_fast):
    re = fmul(a, p, prec, rnd)
    im = fmul(b, p, prec, rnd)
    return re, im

def mpc_mul_int((a, b), n, prec, rnd=round_fast):
    re = fmuli(a, n, prec, rnd)
    im = fmuli(b, n, prec, rnd)
    return re, im

def mpc_div((a, b), (c, d), prec, rnd=round_fast):
    wp = prec + 10
    # mag = c*c + d*d
    mag = fadd(fmul(c, c), fmul(d, d), wp)
    # (a*c+b*d)/mag, (b*c-a*d)/mag
    t = fadd(fmul(a,c), fmul(b,d), wp)
    u = fsub(fmul(b,c), fmul(a,d), wp)
    return fdiv(t,mag,prec,rnd), fdiv(u,mag,prec,rnd)

def mpc_div_mpf((a, b), p, prec, rnd=round_fast):
    re = fdiv(a, p, prec, rnd)
    im = fdiv(b, p, prec, rnd)
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
        re = from_man_exp(re, n*aexp, prec, rnd)
        im = from_man_exp(im, n*bexp, prec, rnd)
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
            im = fsqrt(fneg(a), prec, rnd)
            return (fzero, im)
        else:
            re = fsqrt(a, prec, rnd)
            return (re, fzero)
    wp = prec+20
    if not a[0]:                               # case a positive
        t  = fadd(mpc_abs((a, b), wp), a, wp)  # t = abs(a+bi) + a
        u = fshift(t, -1)                      # u = t/2
        re = fsqrt(u, prec, rnd)               # re = sqrt(u)
        v = fshift(t, 1)                       # v = 2*t
        w  = fsqrt(v, wp)                      # w = sqrt(v)
        im = fdiv(b, w, prec, rnd)             # im = b / w
    else:                                      # case a negative
        t = fsub(mpc_abs((a, b), wp), a, wp)   # t = abs(a+bi) - a
        u = fshift(t, -1)                      # u = t/2
        im = fsqrt(u, prec, rnd)               # im = sqrt(u)
        v = fshift(t, 1)                       # v = 2*t
        w  = fsqrt(v, wp)                      # w = sqrt(v)
        re = fdiv(b, w, prec, rnd)             # re = b/w
        if b[0]:
            re = fneg(re)
            im = fneg(im)
    return re, im

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
    mag = fexp(a, prec+4, rnd)
    c, s = cos_sin(b, prec+4, rnd)
    re = fmul(mag, c, prec, rnd)
    im = fmul(mag, s, prec, rnd)
    return re, im

def mpc_log(z, prec, rnd=round_fast):
    return flog(mpc_abs(z, prec, rnd), prec, rnd), mpc_arg(z, prec, rnd)

def mpc_cos((a, b), prec, rnd=round_fast):
    """Complex cosine. The formula used is cos(a+bi) = cos(a)*cosh(b) -
    sin(a)*sinh(b)*i.

    The same comments apply as for the complex exp: only real
    multiplications are pewrormed, so no cancellation errors are
    possible. The formula is also efficient since we can compute both
    pairs (cos, sin) and (cosh, sinh) in single stwps."""
    if a == fzero:
        return fcosh(b, prec, rnd), fzero
    wp = prec + 6
    c, s = cos_sin(a, wp)
    ch, sh = cosh_sinh(b, wp)
    re = fmul(c, ch, prec, rnd)
    im = fmul(s, sh, prec, rnd)
    return re, fneg(im)

def mpc_sin((a, b), prec, rnd=round_fast):
    """Complex sine. We have sin(a+bi) = sin(a)*cosh(b) +
    cos(a)*sinh(b)*i. See the docstring for mpc_cos for additional
    comments."""
    if a == fzero:
        return fzero, fsinh(b, prec, rnd)
    wp = prec + 6
    c, s = cos_sin(a, wp)
    ch, sh = cosh_sinh(b, wp)
    re = fmul(s, ch, prec, rnd)
    im = fmul(c, sh, prec, rnd)
    return re, im

def mpc_tan(z, prec, rnd=round_fast):
    """Complex tangent. Computed as tan(a+bi) = sin(2a)/M + sinh(2b)/M*i
    where M = cos(2a) + cosh(2b)."""
    a, b = z
    asign, aman, aexp, abc = a
    bsign, bman, bexp, bbc = b
    if b == fzero: return ftan(a, prec, rnd), fzero
    if a == fzero: return fzero, ftanh(b, prec, rnd)
    wp = prec + 15
    a = fshift(a, 1)
    b = fshift(b, 1)
    c, s = cos_sin(a, wp)
    ch, sh = cosh_sinh(b, wp)
    # TODO: handle cancellation when c ~=  -1 and ch ~= 1
    mag = fadd(c, ch, wp)
    re = fdiv(s, mag, prec, rnd)
    im = fdiv(sh, mag, prec, rnd)
    return re, im

def mpc_cosh((a, b), prec, rnd=round_fast):
    """Complex hyperbolic cosine. Computed as cosh(z) = cos(z*i)."""
    return mpc_cos((b, fneg(a)), prec, rnd)

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
    x = fadd(fone, b, wp), fneg(a)
    y = fsub(fone, b, wp), a
    l1 = mpc_log(x, wp)
    l2 = mpc_log(y, wp)
    a, b = mpc_sub(l1, l2, prec, rnd)
    # (I/2) * (a+b*I) = (-b/2 + a/2*I)
    return fneg(fshift(b,-1)), fshift(a,-1)

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
        am = fsub(fone, fabs(a), wp)
        # case abs(a) <= 1 
        if not am[0]:
            if n == 0:
                return facos(a, prec, rnd), fzero
            else:
                return fasin(a, prec, rnd), fzero
        # cases abs(a) > 1
        else:
            # case a < -1
            if a[0]:
                pi = fpi(prec, rnd)
                c = facosh(fneg(a), prec, rnd)
                if n == 0:
                    return pi, fneg(c)
                else:
                    return fneg(fshift(pi, -1)), c
            # case a > 1
            else:
                c = facosh(a, prec, rnd)
                if n == 0:
                    return fzero, c
                else:
                    pi = fpi(prec, rnd)
                    return fshift(pi, -1), fneg(c)
    asign = bsign = 0
    if a[0]:
        a = fneg(a)
        asign = 1
    if b[0]:
        b = fneg(b)
        bsign = 1
    am = fsub(fone, a, wp)
    ap = fadd(fone, a, wp)
    r = fhypot(ap, b, wp)
    s = fhypot(am, b, wp)
    alpha = fshift(fadd(r, s, wp), -1)
    beta = fdiv(a, alpha, wp)
    b2 = fmul(b,b, wp)
    # case beta <= beta_crossover
    if not fsub(beta_crossover, beta, wp)[0]:
        if n == 0:
            re = facos(beta, wp)
        else:
            re = fasin(beta, wp)
    else:
        # to compute the real part in this region use the identity
        # asin(beta) = atan(beta/sqrt(1-beta**2))
        # beta/sqrt(1-beta**2) = (alpha + a) * (alpha - a)
        # alpha + a is numerically accurate; alpha - a can have 
        # cancellations leading to numerical inaccuracies, so rewrite
        # it in differente ways according to the region
        Ax = fadd(alpha, a, wp)
        # case a <= 1
        if not am[0]:
            # c = b*b/(r + (a+1)); d = (s + (1-a))
            # alpha - a = (1/2)*(c + d)
            # case n=0: re = atan(sqrt((1/2) * Ax * (c + d))/a)
            # case n=1: re = atan(a/sqrt((1/2) * Ax * (c + d)))
            c = fdiv(b2, fadd(r, ap, wp), wp)
            d = fadd(s, am, wp)
            re = fshift(fmul(Ax, fadd(c, d, wp), wp), -1)
            if n == 0:
                re = fatan(fdiv(fsqrt(re, wp), a, wp), wp)
            else:
                re = fatan(fdiv(a, fsqrt(re, wp), wp), wp)
        else:
            # c = Ax/(r + (a+1)); d = Ax/(s - (1-a))
            # alpha - a = (1/2)*(c + d)
            # case n = 0: re = atan(b*sqrt(c + d)/2/a)
            # case n = 1: re = atan(a/(b*sqrt(c + d)/2)
            c = fdiv(Ax, fadd(r, ap, wp), wp)
            d = fdiv(Ax, fsub(s, am, wp), wp)
            re = fshift(fadd(c, d, wp), -1)
            re = fmul(b, fsqrt(re, wp), wp)
            if n == 0:
                re = fatan(fdiv(re, a, wp), wp)
            else:
                re = fatan(fdiv(a, re, wp), wp)
    # to compute alpha + sqrt(alpha**2 - 1), if alpha <= alpha_crossover
    # replace it with 1 + Am1 + sqrt(Am1*(alpha+1)))
    # where Am1 = alpha -1
    # if alpha <= alpha_crossover:
    if not fsub(alpha_crossover, alpha, wp)[0]:
        c1 = fdiv(b2, fadd(r, ap, wp), wp)
        # case a < 1
        if fneg(am)[0]:
            # Am1 = (1/2) * (b*b/(r + (a+1)) + b*b/(s + (1-a))
            c2 = fadd(s, am, wp)
            c2 = fdiv(b2, c2, wp)
            Am1 = fshift(fadd(c1, c2, wp), -1)
        else:
            # Am1 = (1/2) * (b*b/(r + (a+1)) + (s - (1-a)))
            c2 = fsub(s, am, wp)
            Am1 = fshift(fadd(c1, c2, wp), -1)
        # im = log(1 + Am1 + sqrt(Am1*(alpha+1)))
        im = fmul(Am1, fadd(alpha, fone, wp), wp)
        im = flog(fadd(fone, fadd(Am1, fsqrt(im, wp), wp), wp), wp)
    else:
        # im = log(alpha + sqrt(alpha*alpha - 1))
        im = fsqrt(fsub(fmul(alpha, alpha, wp), fone, wp), wp)
        im = flog(fadd(alpha, im, wp), wp)
    if asign:
        if n == 0:
            re = fsub(fpi(wp), re, wp)
        else:
            re = fneg(re)
    if not bsign and n == 0:
        im = fneg(im)
    if bsign and n == 1:
        im = fneg(im)
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
    a, b =  mpc_asin((b, fneg(a)), prec, rnd)
    return fneg(b), a

def mpc_acosh(z, prec, rnd=round_fast):
    # acosh(z) = -I * acos(z)   for Im(acos(z)) <= 0
    #            +I * acos(z)   otherwise
    a, b = mpc_acos(z, prec, rnd)
    if b[0] or b == fzero:
        return fneg(b), a
    else:
        return b, fneg(a)

def mpc_atanh(z, prec, rnd=round_fast):
    # atanh(z) = (log(1+z)-log(1-z))/2
    wp = prec + 15
    a = mpc_add(z, mpc_one, wp)
    b = mpc_sub(mpc_one, z, wp)
    a = mpc_log(a, wp)
    b = mpc_log(b, wp)
    return mpc_shift(mpc_sub(a, b, wp), -1)
