"""
Computational functions for interval arithmetic.

"""

from libmpf import (
    ComplexResult,
    round_down, round_up, round_floor, round_ceiling, round_nearest,
    prec_to_dps, repr_dps,
    fnan, finf, fninf, fzero, fhalf, fone, fnone,
    mpf_sign, mpf_lt, mpf_le, mpf_gt, mpf_ge, mpf_eq, mpf_cmp,
    mpf_floor, from_int, to_int, to_str,
    mpf_abs, mpf_neg, mpf_pos, mpf_add, mpf_sub, mpf_mul,
    mpf_div, mpf_shift, mpf_pow_int)

from libelefun import (
    mpf_log, mpf_exp, mpf_sqrt, reduce_angle, calc_cos_sin
)

def mpi_str(s, prec):
    sa, sb = s
    dps = prec_to_dps(prec) + 5
    return "[%s, %s]" % (to_str(sa, dps), to_str(sb, dps))

    #dps = prec_to_dps(prec)
    #m = mpi_mid(s, prec)
    #d = mpf_shift(mpi_delta(s, 20), -1)
    #return "%s +/- %s" % (to_str(m, dps), to_str(d, 3))

def mpi_add(s, t, prec):
    sa, sb = s
    ta, tb = t
    a = mpf_add(sa, ta, prec, round_floor)
    b = mpf_add(sb, tb, prec, round_ceiling)
    if a == fnan: a = fninf
    if b == fnan: b = finf
    return a, b

def mpi_sub(s, t, prec):
    sa, sb = s
    ta, tb = t
    a = mpf_sub(sa, tb, prec, round_floor)
    b = mpf_sub(sb, ta, prec, round_ceiling)
    if a == fnan: a = fninf
    if b == fnan: b = finf
    return a, b

def mpi_delta(s, prec):
    sa, sb = s
    return mpf_sub(sb, sa, prec, round_up)

def mpi_mid(s, prec):
    sa, sb = s
    return mpf_shift(mpf_add(sa, sb, prec, round_nearest), -1)

def mpi_pos(s, prec):
    sa, sb = s
    a = mpf_pos(sa, prec, round_floor)
    b = mpf_pos(sb, prec, round_ceiling)
    return a, b

def mpi_neg(s, prec=None):
    sa, sb = s
    a = mpf_neg(sb, prec, round_floor)
    b = mpf_neg(sa, prec, round_ceiling)
    return a, b

def mpi_abs(s, prec):
    sa, sb = s
    sas = mpf_sign(sa)
    sbs = mpf_sign(sb)
    # Both points nonnegative?
    if sas >= 0:
        a = mpf_pos(sa, prec, round_floor)
        b = mpf_pos(sb, prec, round_ceiling)
    # Upper point nonnegative?
    elif sbs >= 0:
        a = fzero
        negsa = mpf_neg(sa)
        if mpf_lt(negsa, sb):
            b = mpf_pos(sb, prec, round_ceiling)
        else:
            b = mpf_pos(negsa, prec, round_ceiling)
    # Both negative?
    else:
        a = mpf_neg(sb, prec, round_floor)
        b = mpf_neg(sa, prec, round_ceiling)
    return a, b

def mpi_mul(s, t, prec):
    sa, sb = s
    ta, tb = t
    sas = mpf_sign(sa)
    sbs = mpf_sign(sb)
    tas = mpf_sign(ta)
    tbs = mpf_sign(tb)
    if sas == sbs == 0:
        # Should maybe be undefined
        if ta == fninf or tb == finf:
            return fninf, finf
        return fzero, fzero
    if tas == tbs == 0:
        # Should maybe be undefined
        if sa == fninf or sb == finf:
            return fninf, finf
        return fzero, fzero
    if sas >= 0:
        # positive * positive
        if tas >= 0:
            a = mpf_mul(sa, ta, prec, round_floor)
            b = mpf_mul(sb, tb, prec, round_ceiling)
            if a == fnan: a = fzero
            if b == fnan: b = finf
        # positive * negative
        elif tbs <= 0:
            a = mpf_mul(sb, ta, prec, round_floor)
            b = mpf_mul(sa, tb, prec, round_ceiling)
            if a == fnan: a = fninf
            if b == fnan: b = fzero
        # positive * both signs
        else:
            a = mpf_mul(sb, ta, prec, round_floor)
            b = mpf_mul(sb, tb, prec, round_ceiling)
            if a == fnan: a = fninf
            if b == fnan: b = finf
    elif sbs <= 0:
        # negative * positive
        if tas >= 0:
            a = mpf_mul(sa, tb, prec, round_floor)
            b = mpf_mul(sb, ta, prec, round_ceiling)
            if a == fnan: a = fninf
            if b == fnan: b = fzero
        # negative * negative
        elif tbs <= 0:
            a = mpf_mul(sb, tb, prec, round_floor)
            b = mpf_mul(sa, ta, prec, round_ceiling)
            if a == fnan: a = fzero
            if b == fnan: b = finf
        # negative * both signs
        else:
            a = mpf_mul(sa, tb, prec, round_floor)
            b = mpf_mul(sa, ta, prec, round_ceiling)
            if a == fnan: a = fninf
            if b == fnan: b = finf
    else:
        # General case: perform all cross-multiplications and compare
        # Since the multiplications can be done exactly, we need only
        # do 4 (instead of 8: two for each rounding mode)
        cases = [mpf_mul(sa, ta), mpf_mul(sa, tb), mpf_mul(sb, ta), mpf_mul(sb, tb)]
        if fnan in cases:
            a, b = (fninf, finf)
        else:
            cases = sorted(cases, cmp=mpf_cmp)
            a = mpf_pos(cases[0], prec, round_floor)
            b = mpf_pos(cases[-1], prec, round_ceiling)
    return a, b

def mpi_div(s, t, prec):
    sa, sb = s
    ta, tb = t
    sas = mpf_sign(sa)
    sbs = mpf_sign(sb)
    tas = mpf_sign(ta)
    tbs = mpf_sign(tb)
    # 0 / X
    if sas == sbs == 0:
        # 0 / <interval containing 0>
        if (tas < 0 and tbs > 0) or (tas == 0 or tbs == 0):
            return fninf, finf
        return fzero, fzero
    # Denominator contains both negative and positive numbers;
    # this should properly be a multi-interval, but the closest
    # match is the entire (extended) real line
    if tas < 0 and tbs > 0:
        return fninf, finf
    # Assume denominator to be nonnegative
    if tas < 0:
        return mpi_div(mpi_neg(s), mpi_neg(t), prec)
    # Division by zero
    # XXX: make sure all results make sense
    if tas == 0:
        # Numerator contains both signs?
        if sas < 0 and sbs > 0:
            return fninf, finf
        if tas == tbs:
            return fninf, finf
        # Numerator positive?
        if sas >= 0:
            a = mpf_div(sa, tb, prec, round_floor)
            b = finf
        if sbs <= 0:
            a = fninf
            b = mpf_div(sb, tb, prec, round_ceiling)
    # Division with positive denominator
    # We still have to handle nans resulting from inf/0 or inf/inf
    else:
        # Nonnegative numerator
        if sas >= 0:
            a = mpf_div(sa, tb, prec, round_floor)
            b = mpf_div(sb, ta, prec, round_ceiling)
            if a == fnan: a = fzero
            if b == fnan: b = finf
        # Nonpositive numerator
        elif sbs <= 0:
            a = mpf_div(sa, ta, prec, round_floor)
            b = mpf_div(sb, tb, prec, round_ceiling)
            if a == fnan: a = fninf
            if b == fnan: b = fzero
        # Numerator contains both signs?
        else:
            a = mpf_div(sa, ta, prec, round_floor)
            b = mpf_div(sb, ta, prec, round_ceiling)
            if a == fnan: a = fninf
            if b == fnan: b = finf
    return a, b

def mpi_exp(s, prec):
    sa, sb = s
    # exp is monotonous
    a = mpf_exp(sa, prec, round_floor)
    b = mpf_exp(sb, prec, round_ceiling)
    return a, b

def mpi_log(s, prec):
    sa, sb = s
    # log is monotonous
    a = mpf_log(sa, prec, round_floor)
    b = mpf_log(sb, prec, round_ceiling)
    return a, b

def mpi_sqrt(s, prec):
    sa, sb = s
    # sqrt is monotonous
    a = mpf_sqrt(sa, prec, round_floor)
    b = mpf_sqrt(sb, prec, round_ceiling)
    return a, b

def mpi_pow_int(s, n, prec):
    sa, sb = s
    if n < 0:
        return mpi_div((fone, fone), mpi_pow_int(s, -n, prec+20), prec)
    if n == 0:
        return (fone, fone)
    if n == 1:
        return s
    # Odd -- signs are preserved
    if n & 1:
        a = mpf_pow_int(sa, n, prec, round_floor)
        b = mpf_pow_int(sb, n, prec, round_ceiling)
    # Even -- important to ensure positivity
    else:
        sas = mpf_sign(sa)
        sbs = mpf_sign(sb)
        # Nonnegative?
        if sas >= 0:
            a = mpf_pow_int(sa, n, prec, round_floor)
            b = mpf_pow_int(sb, n, prec, round_ceiling)
        # Nonpositive?
        elif sbs <= 0:
            a = mpf_pow_int(sb, n, prec, round_floor)
            b = mpf_pow_int(sa, n, prec, round_ceiling)
        # Mixed signs?
        else:
            a = fzero
            # max(-a,b)**n
            sa = mpf_neg(sa)
            if mpf_ge(sa, sb):
                b = mpf_pow_int(sa, n, prec, round_ceiling)
            else:
                b = mpf_pow_int(sb, n, prec, round_ceiling)
    return a, b

def mpi_pow(s, t, prec):
    ta, tb = t
    if ta == tb and ta not in (finf, fninf):
        if ta == from_int(to_int(ta)):
            return mpi_pow_int(s, to_int(ta), prec)
        if ta == fhalf:
            return mpi_sqrt(s, prec)
    u = mpi_log(s, prec + 20)
    v = mpi_mul(u, t, prec + 20)
    return mpi_exp(v, prec)

def MIN(x, y):
    if mpf_le(x, y):
        return x
    return y

def MAX(x, y):
    if mpf_ge(x, y):
        return x
    return y

def mpi_cos_sin(x, prec):
    a, b = x
    # Guaranteed to contain both -1 and 1
    if finf in (a, b) or fninf in (a, b):
        return (fnone, fone), (fnone, fone)
    y, yswaps, yn = reduce_angle(a, prec+20)
    z, zswaps, zn = reduce_angle(b, prec+20)
    # Guaranteed to contain both -1 and 1
    if zn - yn >= 4:
        return (fnone, fone), (fnone, fone)
    # Both points in the same quadrant -- cos and sin both strictly monotonous
    if yn == zn:
        m = yn % 4
        if m == 0:
            cb, sa = calc_cos_sin(0, y, yswaps, prec, round_ceiling, round_floor)
            ca, sb = calc_cos_sin(0, z, zswaps, prec, round_floor, round_ceiling)
        if m == 1:
            cb, sb = calc_cos_sin(0, y, yswaps, prec, round_ceiling, round_ceiling)
            ca, sa = calc_cos_sin(0, z, zswaps, prec, round_floor, round_ceiling)
        if m == 2:
            ca, sb = calc_cos_sin(0, y, yswaps, prec, round_floor, round_ceiling)
            cb, sa = calc_cos_sin(0, z, zswaps, prec, round_ceiling, round_floor)
        if m == 3:
            ca, sa = calc_cos_sin(0, y, yswaps, prec, round_floor, round_floor)
            cb, sb = calc_cos_sin(0, z, zswaps, prec, round_ceiling, round_ceiling)
        return (ca, cb), (sa, sb)
    # Intervals spanning multiple quadrants
    yn %= 4
    zn %= 4
    case = (yn, zn)
    if case == (0, 1):
        cb, sy = calc_cos_sin(0, y, yswaps, prec, round_ceiling, round_floor)
        ca, sz = calc_cos_sin(0, z, zswaps, prec, round_floor, round_floor)
        return (ca, cb), (MIN(sy, sz), fone)
    if case == (3, 0):
        cy, sa = calc_cos_sin(0, y, yswaps, prec, round_floor, round_floor)
        cz, sb = calc_cos_sin(0, z, zswaps, prec, round_floor, round_ceiling)
        return (MIN(cy, cz), fone), (sa, sb)


    raise NotImplementedError("cos/sin spanning multiple quadrants")

def mpi_cos(x, prec):
    return mpi_cos_sin(x, prec)[0]

def mpi_sin(x, prec):
    return mpi_cos_sin(x, prec)[1]

def mpi_tan(x, prec):
    cos, sin = mpi_cos_sin(x, prec+20)
    return mpi_div(sin, cos, prec)

def mpi_cot(x, prec):
    cos, sin = mpi_cos_sin(x, prec+20)
    return mpi_div(cos, sin, prec)
