#!/usr/bin/env python
"""
    elliptic.py

    Implements the Jacobi theta and Jacobi elliptic functions, using
    arbitrary precision math library

    Author of the first version: M.T. Taschuk

    References:

    [1] Abramowitz & Stegun. 'Handbook of Mathematical Functions, 9th Ed.',
        (Dover duplicate of 1972 edition)
    [2] Whittaker 'A Course of Modern Analysis, 4th Ed.', 1946,
        Cambridge Univeristy Press

"""
from mptypes import (mpf, mpc, mp, mpmathify, eps, one, zero, j)
from functions import (pi, sqrt, cos, sin, exp, log, tanh, ellipk,
                       sech, nthroot)
from libmpf import to_fixed, MP_ZERO, mpf_shift, from_man_exp
from libelefun import cos_sin

# The series for the Jacobi theta functions converge for |q| < 1;
# in the current implementation they throw a ValueError for
# abs(q) > Q_LIM
Q_LIM = 1 - 10**-7

def calculate_nome(k):
    """
    Calculate the nome, q, from the value for k.

    Useful factoids:

    k**2 = m;   m is used in Abramowitz
    """
    k = mpmathify(k)

    if abs(k) > one:             # range error
        raise ValueError

    if k == zero:
        return zero
    elif k == one:
        return one
    else:
        kprimesquared = one - k**2
        kprime = sqrt(kprimesquared)
        top = ellipk(kprimesquared)
        bottom = ellipk(k**2)

        argument = -pi*top/bottom

        nome = exp(argument)
        return nome

# XXX: unused
def calculate_k(q):
    """
    Calculates the value of k for a particular nome, q,
    using Jacobi theta functions.
    """

    q = mpmathify(q)

    v2 = jtheta(2, 0, q)
    v3 = jtheta(3, 0, q)
    m = v2**2/v3**2
    return m


def _jacobi_theta2(z, q):
    extra1 = 10
    extra2 = 20
    # the loops below break when the fixed precision quantities
    # a and b go to zero;
    # right shifting small negative numbers by wp one obtains -1, not zero,
    # so the condition a**2 + b**2 > MIN is used to break the loops.
    MIN = 2
    if z == zero:
        if isinstance(q, mpf):
            wp = mp.prec + extra1
            x = to_fixed(q._mpf_, wp)
            x2 = (x*x) >> wp
            a = b = x2
            s = x2
            while abs(a) > MIN:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                s += a
            s = (1 << (wp+1)) + (s << 1)
            s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
        else:
            wp = mp.prec + extra1
            xre, xim = q._mpc_
            xre = to_fixed(xre, wp)
            xim = to_fixed(xim, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = x2re
            aim = bim = x2im
            sre = (1<<wp) + are
            sim = aim
            while are**2 + aim**2 > MIN:
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                sre += are
                sim += aim
            sre = (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
    else:
        if isinstance(q, mpf) and isinstance(z, mpf):
            wp = mp.prec + extra1
            x = to_fixed(q._mpf_, wp)
            x2 = (x*x) >> wp
            a = b = x2
            c1, s1 = cos_sin(z._mpf_, wp)
            cn = c1 = to_fixed(c1, wp)
            sn = s1 = to_fixed(s1, wp)
            c2 = (c1*c1 - s1*s1) >> wp
            s2 = (c1 * s1) >> (wp - 1)
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
            s = c1 + ((a * cn) >> wp)
            while abs(a) > MIN:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
                s += (a * cn) >> wp
            s = (s << 1)
            s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
            s *= nthroot(q, 4)
            return s
        # case z real, q complex
        elif isinstance(z, mpf):
            wp = mp.prec + extra2
            xre, xim = q._mpc_
            xre = to_fixed(xre, wp)
            xim = to_fixed(xim, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = x2re
            aim = bim = x2im
            c1, s1 = cos_sin(z._mpf_, wp)
            cn = c1 = to_fixed(c1, wp)
            sn = s1 = to_fixed(s1, wp)
            c2 = (c1*c1 - s1*s1) >> wp
            s2 = (c1 * s1) >> (wp - 1)
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
            sre = c1 + ((are * cn) >> wp)
            sim = ((aim * cn) >> wp)
            while are**2 + aim**2 > MIN:
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp

                sre += ((are * cn) >> wp)
                sim += ((aim * cn) >> wp)
            sre = (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
        #case z complex, q real
        elif isinstance(q, mpf):
            wp = mp.prec + extra2
            x = to_fixed(q._mpf_, wp)
            x2 = (x*x) >> wp
            a = b = x2
            prec0 = mp.prec
            mp.prec = wp
            c1 = cos(z)
            s1 = sin(z)
            mp.prec = prec0
            cnre = c1re = to_fixed(c1.real._mpf_, wp)
            cnim = c1im = to_fixed(c1.imag._mpf_, wp)
            snre = s1re = to_fixed(s1.real._mpf_, wp)
            snim = s1im = to_fixed(s1.imag._mpf_, wp)
            #c2 = (c1*c1 - s1*s1) >> wp
            c2re = (c1re*c1re - c1im*c1im - s1re*s1re + s1im*s1im) >> wp
            c2im = (c1re*c1im - s1re*s1im) >> (wp - 1)
            #s2 = (c1 * s1) >> (wp - 1)
            s2re = (c1re*s1re - c1im*s1im) >> (wp - 1)
            s2im = (c1re*s1im + c1im*s1re) >> (wp - 1)
            #cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
            t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
            t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
            t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
            t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4

            sre = c1re + ((a * cnre) >> wp)
            sim = c1im + ((a * cnim) >> wp)
            while abs(a) > MIN:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
                t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
                t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
                t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
                cnre = t1
                cnim = t2
                snre = t3
                snim = t4
                sre += ((a * cnre) >> wp)
                sim += ((a * cnim) >> wp)
            sre = (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
        # case z and q complex
        else:
            wp = mp.prec + extra2
            xre, xim = q._mpc_
            xre = to_fixed(xre, wp)
            xim = to_fixed(xim, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = x2re
            aim = bim = x2im
            prec0 = mp.prec
            mp.prec = wp
            # cos(z), siz(z) with z complex
            c1 = cos(z)
            s1 = sin(z)
            mp.prec = prec0
            cnre = c1re = to_fixed(c1.real._mpf_, wp)
            cnim = c1im = to_fixed(c1.imag._mpf_, wp)
            snre = s1re = to_fixed(s1.real._mpf_, wp)
            snim = s1im = to_fixed(s1.imag._mpf_, wp)
            c2re = (c1re*c1re - c1im*c1im - s1re*s1re + s1im*s1im) >> wp
            c2im = (c1re*c1im - s1re*s1im) >> (wp - 1)
            s2re = (c1re*s1re - c1im*s1im) >> (wp - 1)
            s2im = (c1re*s1im + c1im*s1re) >> (wp - 1)
            t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
            t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
            t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
            t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            n = 1
            termre = c1re
            termim = c1im
            sre = c1re + ((are * cnre - aim * cnim) >> wp)
            sim = c1im + ((are * cnim + aim * cnre) >> wp)

            n = 3
            termre = ((are * cnre - aim * cnim) >> wp)
            termim = ((are * cnim + aim * cnre) >> wp)
            sre = c1re + ((are * cnre - aim * cnim) >> wp)
            sim = c1im + ((are * cnim + aim * cnre) >> wp)

            n = 5
            while are**2 + aim**2 > MIN:
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                #cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
                t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
                t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
                t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
                t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
                cnre = t1
                cnim = t2
                snre = t3
                snim = t4
                termre = ((are * cnre - aim * cnim) >> wp)
                termim = ((aim * cnre + are * cnim) >> wp)
                sre += ((are * cnre - aim * cnim) >> wp)
                sim += ((aim * cnre + are * cnim) >> wp)
                n += 2
            sre = (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
    s *= nthroot(q, 4)
    return s

def _djacobi_theta2(z, q, nd):
    MIN = 2
    extra1 = 10
    extra2 = 20
    if isinstance(q, mpf) and isinstance(z, mpf):
        wp = mp.prec + extra1
        x = to_fixed(q._mpf_, wp)
        x2 = (x*x) >> wp
        a = b = x2
        c1, s1 = cos_sin(z._mpf_, wp)
        cn = c1 = to_fixed(c1, wp)
        sn = s1 = to_fixed(s1, wp)
        c2 = (c1*c1 - s1*s1) >> wp
        s2 = (c1 * s1) >> (wp - 1)
        cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
        if (nd&1):
            s = s1 + ((a * sn * 3**nd) >> wp)
        else:
            s = c1 + ((a * cn * 3**nd) >> wp)
        n = 2
        while abs(a) > MIN:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
            if nd&1:
                s += (a * sn * (2*n+1)**nd) >> wp
            else:
                s += (a * cn * (2*n+1)**nd) >> wp
            n += 1
        s = -(s << 1)
        s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
        # case z real, q complex
    elif isinstance(z, mpf):
        wp = mp.prec + extra2
        xre, xim = q._mpc_
        xre = to_fixed(xre, wp)
        xim = to_fixed(xim, wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = x2re
        aim = bim = x2im
        c1, s1 = cos_sin(z._mpf_, wp)
        cn = c1 = to_fixed(c1, wp)
        sn = s1 = to_fixed(s1, wp)
        c2 = (c1*c1 - s1*s1) >> wp
        s2 = (c1 * s1) >> (wp - 1)
        cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
        if (nd&1):
            sre = s1 + ((are * sn * 3**nd) >> wp)
            sim = ((aim * sn * 3**nd) >> wp)
        else:
            sre = c1 + ((are * cn * 3**nd) >> wp)
            sim = ((aim * cn * 3**nd) >> wp)
        n = 5
        while are**2 + aim**2 > MIN:
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp

            if (nd&1):
                sre += ((are * sn * n**nd) >> wp)
                sim += ((aim * sn * n**nd) >> wp)
            else:
                sre += ((are * cn * n**nd) >> wp)
                sim += ((aim * cn * n**nd) >> wp)
            n += 2
        sre = -(sre << 1)
        sim = -(sim << 1)
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
    #case z complex, q real
    elif isinstance(q, mpf):
        wp = mp.prec + extra2
        x = to_fixed(q._mpf_, wp)
        x2 = (x*x) >> wp
        a = b = x2
        prec0 = mp.prec
        mp.prec = wp
        c1 = cos(z)
        s1 = sin(z)
        mp.prec = prec0
        cnre = c1re = to_fixed(c1.real._mpf_, wp)
        cnim = c1im = to_fixed(c1.imag._mpf_, wp)
        snre = s1re = to_fixed(s1.real._mpf_, wp)
        snim = s1im = to_fixed(s1.imag._mpf_, wp)
        #c2 = (c1*c1 - s1*s1) >> wp
        c2re = (c1re*c1re - c1im*c1im - s1re*s1re + s1im*s1im) >> wp
        c2im = (c1re*c1im - s1re*s1im) >> (wp - 1)
        #s2 = (c1 * s1) >> (wp - 1)
        s2re = (c1re*s1re - c1im*s1im) >> (wp - 1)
        s2im = (c1re*s1im + c1im*s1re) >> (wp - 1)
        #cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
        t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
        t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
        t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
        t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
        cnre = t1
        cnim = t2
        snre = t3
        snim = t4

        if (nd&1):
            sre = s1re + ((a * snre * 3**nd) >> wp)
            sim = s1im + ((a * snim * 3**nd) >> wp)
        else:
            sre = c1re + ((a * cnre * 3**nd) >> wp)
            sim = c1im + ((a * cnim * 3**nd) >> wp)
        n = 5
        while abs(a) > MIN:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
            t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
            t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
            t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            if (nd&1):
                sre += ((a * snre * n**nd) >> wp)
                sim += ((a * snim * n**nd) >> wp)
            else:
                sre += ((a * cnre * n**nd) >> wp)
                sim += ((a * cnim * n**nd) >> wp)
            n += 2
        sre = -(sre << 1)
        sim = -(sim << 1)
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
    # case z and q complex
    else:
        wp = mp.prec + extra2
        xre, xim = q._mpc_
        xre = to_fixed(xre, wp)
        xim = to_fixed(xim, wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = x2re
        aim = bim = x2im
        prec0 = mp.prec
        mp.prec = wp
        # cos(2*z), siz(2*z) with z complex
        c1 = cos(z)
        s1 = sin(z)
        mp.prec = prec0
        cnre = c1re = to_fixed(c1.real._mpf_, wp)
        cnim = c1im = to_fixed(c1.imag._mpf_, wp)
        snre = s1re = to_fixed(s1.real._mpf_, wp)
        snim = s1im = to_fixed(s1.imag._mpf_, wp)
        c2re = (c1re*c1re - c1im*c1im - s1re*s1re + s1im*s1im) >> wp
        c2im = (c1re*c1im - s1re*s1im) >> (wp - 1)
        s2re = (c1re*s1re - c1im*s1im) >> (wp - 1)
        s2im = (c1re*s1im + c1im*s1re) >> (wp - 1)
        t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
        t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
        t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
        t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
        cnre = t1
        cnim = t2
        snre = t3
        snim = t4
        if (nd&1):
            sre = s1re + (((are * snre - aim * snim) * 3**nd) >> wp)
            sim = s1im + (((are * snim + aim * snre)* 3**nd) >> wp)
        else:
            sre = c1re + (((are * cnre - aim * cnim) * 3**nd) >> wp)
            sim = c1im + (((are * cnim + aim * cnre)* 3**nd) >> wp)
        n = 5
        while are**2 + aim**2 > MIN:
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            #cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
            t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
            t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
            t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
            t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            if (nd&1):
                sre += (((are * snre - aim * snim) * n**nd) >> wp)
                sim += (((aim * snre + are * snim) * n**nd) >> wp)
            else:
                sre += (((are * cnre - aim * cnim) * n**nd) >> wp)
                sim += (((aim * cnre + are * cnim) * n**nd) >> wp)
            n += 2
        sre = -(sre << 1)
        sim = -(sim << 1)
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
    s *= nthroot(q, 4)
    if (nd&1):
        return (-1)**(nd//2) * s
    else:
        return (-1)**(1 + nd//2) * s

def _jacobi_theta3(z, q):
    extra1 = 10
    extra2 = 20
    MIN = 2
    if z == zero:
        if isinstance(q, mpf):
            wp = mp.prec + extra1
            x = to_fixed(q._mpf_, wp)
            s = x
            a = b = x
            x2 = (x*x) >> wp
            while abs(a) > MIN:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                s += a
            s = (1 << wp) + (s << 1)
            s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
            return s
        else:
            wp = mp.prec + extra1
            xre, xim = q._mpc_
            xre = to_fixed(xre, wp)
            xim = to_fixed(xim, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            sre = are = bre = xre
            sim = aim = bim = xim
            while are**2 + aim**2 > MIN:
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                sre += are
                sim += aim
            sre = (1 << wp) + (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
            return s
    else:
        if isinstance(q, mpf) and isinstance(z, mpf):
            s = MP_ZERO
            wp = mp.prec + extra1
            x = to_fixed(q._mpf_, wp)
            a = b = x
            x2 = (x*x) >> wp
            c1, s1 = cos_sin(mpf_shift(z._mpf_, 1), wp)
            c1 = to_fixed(c1, wp)
            s1 = to_fixed(s1, wp)
            cn = c1
            sn = s1
            s += (a * cn) >> wp
            while abs(a) > MIN:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
                s += (a * cn) >> wp
            s = (1 << wp) + (s << 1)
            s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
            return s
        # case z real, q complex
        elif isinstance(z, mpf):
            wp = mp.prec + extra2
            xre, xim = q._mpc_
            xre = to_fixed(xre, wp)
            xim = to_fixed(xim, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = xre
            aim = bim = xim
            c1, s1 = cos_sin(mpf_shift(z._mpf_, 1), wp)
            c1 = to_fixed(c1, wp)
            s1 = to_fixed(s1, wp)
            cn = c1
            sn = s1
            sre = (are * cn) >> wp
            sim = (aim * cn) >> wp
            while are**2 + aim**2 > MIN:
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp

                sre += (are * cn) >> wp
                sim += (aim * cn) >> wp
            sre = (1 << wp) + (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
            return s
        #case z complex, q real
        elif isinstance(q, mpf):
            wp = mp.prec + extra2
            x = to_fixed(q._mpf_, wp)
            a = b = x
            x2 = (x*x) >> wp
            prec0 = mp.prec
            mp.prec = wp
            c1 = cos(2*z)
            s1 = sin(2*z)
            mp.prec = prec0
            cnre = c1re = to_fixed(c1.real._mpf_, wp)
            cnim = c1im = to_fixed(c1.imag._mpf_, wp)
            snre = s1re = to_fixed(s1.real._mpf_, wp)
            snim = s1im = to_fixed(s1.imag._mpf_, wp)
            sre = (a * cnre) >> wp
            sim = (a * cnim) >> wp
            while abs(a) > MIN:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
                t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
                t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
                t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
                cnre = t1
                cnim = t2
                snre = t3
                snim = t4
                sre += (a * cnre) >> wp
                sim += (a * cnim) >> wp
            sre = (1 << wp) + (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
            return s
        # case z and q complex
        else:
            wp = mp.prec + extra2
            xre, xim = q._mpc_
            xre = to_fixed(xre, wp)
            xim = to_fixed(xim, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = xre
            aim = bim = xim
            prec0 = mp.prec
            mp.prec = wp
            # cos(2*z), sin(2*z) with z complex
            c1 = cos(2*z)
            s1 = sin(2*z)
            mp.prec = prec0
            cnre = c1re = to_fixed(c1.real._mpf_, wp)
            cnim = c1im = to_fixed(c1.imag._mpf_, wp)
            snre = s1re = to_fixed(s1.real._mpf_, wp)
            snim = s1im = to_fixed(s1.imag._mpf_, wp)
            sre = (are * cnre - aim * cnim) >> wp
            sim = (aim * cnre + are * cnim) >> wp
            while are**2 + aim**2 > MIN:
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
                t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
                t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
                t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
                cnre = t1
                cnim = t2
                snre = t3
                snim = t4
                sre += (are * cnre - aim * cnim) >> wp
                sim += (aim * cnre + are * cnim) >> wp
            sre = (1 << wp) + (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
            return s

def _djacobi_theta3(z, q, nd):
    """nd=1,2,3 order of the derivative with respect to z"""
    MIN = 2
    extra1 = 10
    extra2 = 20
    if isinstance(q, mpf) and isinstance(z, mpf):
        s = MP_ZERO
        wp = mp.prec + extra1
        x = to_fixed(q._mpf_, wp)
        a = b = x
        x2 = (x*x) >> wp
        c1, s1 = cos_sin(mpf_shift(z._mpf_, 1), wp)
        c1 = to_fixed(c1, wp)
        s1 = to_fixed(s1, wp)
        cn = c1
        sn = s1
        if (nd&1):
            s += (a * sn) >> wp
        else:
            s += (a * cn) >> wp
        n = 2
        while abs(a) > MIN:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
            if nd&1:
                s += (a * sn * n**nd) >> wp
            else:
                s += (a * cn * n**nd) >> wp
            n += 1
        s = -(s << (nd+1))
        s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
    # case z real, q complex
    elif isinstance(z, mpf):
        wp = mp.prec + extra2
        xre, xim = q._mpc_
        xre = to_fixed(xre, wp)
        xim = to_fixed(xim, wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = xre
        aim = bim = xim
        c1, s1 = cos_sin(mpf_shift(z._mpf_, 1), wp)
        c1 = to_fixed(c1, wp)
        s1 = to_fixed(s1, wp)
        cn = c1
        sn = s1
        if (nd&1):
            sre = (are * sn) >> wp
            sim = (aim * sn) >> wp
        else:
            sre = (are * cn) >> wp
            sim = (aim * cn) >> wp
        n = 2
        while are**2 + aim**2 > MIN:
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
            if nd&1:
                sre += (are * sn * n**nd) >> wp
                sim += (aim * sn * n**nd) >> wp
            else:
                sre += (are * cn * n**nd) >> wp
                sim += (aim * cn * n**nd) >> wp
            n += 1
        sre = -(sre << (nd+1))
        sim = -(sim << (nd+1))
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
    #case z complex, q real
    elif isinstance(q, mpf):
        wp = mp.prec + extra2
        x = to_fixed(q._mpf_, wp)
        a = b = x
        x2 = (x*x) >> wp
        prec0 = mp.prec
        mp.prec = wp
        c1 = cos(2*z)
        s1 = sin(2*z)
        mp.prec = prec0
        cnre = c1re = to_fixed(c1.real._mpf_, wp)
        cnim = c1im = to_fixed(c1.imag._mpf_, wp)
        snre = s1re = to_fixed(s1.real._mpf_, wp)
        snim = s1im = to_fixed(s1.imag._mpf_, wp)
        if (nd&1):
            sre = (a * snre) >> wp
            sim = (a * snim) >> wp
        else:
            sre = (a * cnre) >> wp
            sim = (a * cnim) >> wp
        n = 2
        while abs(a) > MIN:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
            t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
            t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
            t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            if (nd&1):
                sre += (a * snre * n**nd) >> wp
                sim += (a * snim * n**nd) >> wp
            else:
                sre += (a * cnre * n**nd) >> wp
                sim += (a * cnim * n**nd) >> wp
            n += 1
        sre = -(sre << (nd+1))
        sim = -(sim << (nd+1))
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
    # case z and q complex
    else:
        wp = mp.prec + extra2
        xre, xim = q._mpc_
        xre = to_fixed(xre, wp)
        xim = to_fixed(xim, wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = xre
        aim = bim = xim
        prec0 = mp.prec
        mp.prec = wp
        # cos(2*z), sin(2*z) with z complex
        c1 = cos(2*z)
        s1 = sin(2*z)
        mp.prec = prec0
        cnre = c1re = to_fixed(c1.real._mpf_, wp)
        cnim = c1im = to_fixed(c1.imag._mpf_, wp)
        snre = s1re = to_fixed(s1.real._mpf_, wp)
        snim = s1im = to_fixed(s1.imag._mpf_, wp)
        if (nd&1):
            sre = (are * snre - aim * snim) >> wp
            sim = (aim * snre + are * snim) >> wp
        else:
            sre = (are * cnre - aim * cnim) >> wp
            sim = (aim * cnre + are * cnim) >> wp
        n = 2
        while are**2 + aim**2 > MIN:
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
            t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
            t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
            t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            if(nd&1):
                sre += ((are * snre - aim * snim) * n**nd) >> wp
                sim += ((aim * snre + are * snim) * n**nd) >> wp
            else:
                sre += ((are * cnre - aim * cnim) * n**nd) >> wp
                sim += ((aim * cnre + are * cnim) * n**nd) >> wp
            n += 1
        sre = -(sre << (nd+1))
        sim = -(sim << (nd+1))
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
    if (nd&1):
        return (-1)**(nd//2) * s
    else:
        return (-1)**(1 + nd//2) * s

def _jacobi_theta2a(z, q):
    """
    case z.imag != 0
    theta(2, z, q) =
    q**1/4 * Sum(q**(n*n + n) * exp(j*(2*n + 1)*z), n=-inf, inf)
    max term for minimum (2*n+1)*log(q).real - 2* z.imag
    n0 = int(z.imag/log(q).real - 1/2)
    theta(2, z, q) =
    q**1/4 * Sum(q**(n*n + n) * exp(j*(2*n + 1)*z), n=n0, inf) +
    q**1/4 * Sum(q**(n*n + n) * exp(j*(2*n + 1)*z), n, n0-1, -inf)
    """
    n = n0 = int(z.imag/log(q).real - 1/2)
    e2 = exp(2*j*z)
    e = e0 = exp(j*(2*n + 1)*z)
    a = q**(n*n + n)
    # leading term
    term = a * e
    s = term
    eps1 = eps*abs(term)
    while 1:
        n += 1
        e = e * e2
        term = q**(n*n + n) * e
        if abs(term) < eps1:
            break
        s += term
    e = e0
    e2 = exp(-2*j*z)
    n = n0
    while 1:
        n -= 1
        e = e * e2
        term = q**(n*n + n) * e
        if abs(term) < eps1:
            break
        s += term
    s = s * nthroot(q, 4)
    return s

def _jacobi_theta3a(z, q):
    """
    case z.imag != 0
    theta3(z, q) = Sum(q**(n*n) * exp(j*2*n*z), n, -inf, inf)
    max term for n*abs(log(q).real) + z.imag ~= 0
    n0 = int(- z.imag/abs(log(q).real))
    """
    n = n0 = int(- z.imag/abs(log(q).real))
    e2 = exp(2*j*z)
    e = e0 = exp(j*2*n*z)
    s = term = q**(n*n) * e
    eps1 = eps*abs(term)
    while 1:
        n += 1
        e = e * e2
        term = q**(n*n) * e
        if abs(term) < eps1:
            break
        s += term
    e = e0
    e2 = exp(-2*j*z)
    n = n0
    while 1:
        n -= 1
        e = e * e2
        term = q**(n*n) * e
        if abs(term) < eps1:
            break
        s += term
    return s

def _djacobi_theta2a(z, q, nd):
    """
    case z.imag != 0
    dtheta(2, z, q, nd) =
    j* q**1/4 * Sum(q**(n*n + n) * (2*n+1)*exp(j*(2*n + 1)*z), n=-inf, inf)
    max term for (2*n0+1)*log(q).real - 2* z.imag ~= 0
    n0 = int(z.imag/log(q).real - 1/2)
    """
    n = n0 = int(z.imag/log(q).real - 1/2)
    e2 = exp(2*j*z)
    e = e0 = exp(j*(2*n + 1)*z)
    a = q**(n*n + n)
    # leading term
    term = (2*n+1)**nd * a * e
    s = term
    eps1 = eps*abs(term)
    while 1:
        n += 1
        e = e * e2
        term = (2*n+1)**nd * q**(n*n + n) * e
        if abs(term) < eps1:
            break
        s += term
    e = e0
    e2 = exp(-2*j*z)
    n = n0
    while 1:
        n -= 1
        e = e * e2
        term = (2*n+1)**nd * q**(n*n + n) * e
        if abs(term) < eps1:
            break
        s += term
    return j**nd * s * nthroot(q, 4)

def _djacobi_theta3a(z, q, nd):
    """
    case z.imag != 0
    djtheta3(z, q, nd) = (2*j)**nd *
      Sum(q**(n*n) * n**nd * exp(j*2*n*z), n, -inf, inf)
    max term for minimum n*abs(log(q).real) + z.imag
    """
    n = n0 = int(-z.imag/abs(log(q).real))
    e2 = exp(2*j*z)
    e = e0 = exp(j*2*n*z)
    a = q**(n*n) * e
    s = term = n**nd * a
    if n != 0:
        eps1 = eps*abs(term)
    else:
        eps1 = eps*abs(a)
    while 1:
        n += 1
        e = e * e2
        a = q**(n*n) * e
        term = n**nd * a
        if n != 0:
            aterm = abs(term)
        else:
            aterm = abs(a)
        if aterm < eps1:
            break
        s += term
    e = e0
    e2 = exp(-2*j*z)
    n = n0
    while 1:
        n -= 1
        e = e * e2
        a = q**(n*n) * e
        term = n**nd * a
        if n != 0:
            aterm = abs(term)
        else:
            aterm = abs(a)
        if aterm < eps1:
            break
        s += term
    return (2*j)**nd * s




def jtheta(n, z, q):
    r"""
    Computes the Jacobi theta function `\vartheta_n(z, q)`, where
    `n = 1, 2, 3, 4`. The theta functions are functions of two
    variables:

    * `z` is the *argument*, an arbitrary real or complex number

    * `q` is the *nome*, which must be a real or complex number
      in the unit disk (i.e. `|q| < 1`)

    One also commonly encounters the notation `\vartheta_n(z, \tau)`
    in the literature. The variable `\tau` is called the *parameter*
    and can be converted to a nome using the formula
    `q = \exp(i \pi \tau)`. Note the condition `|q| < 1` requires
    `\Im(\tau) > 0`; i.e. Jacobi theta functions are defined for
    `\tau` in the upper half plane.

    Other notations are also in use. For example, some authors use
    the single-argument form `\vartheta_n(x)`. Depending on context,
    this can mean ``jtheta(n, 0, x)``, ``jtheta(n, x, q)``, or possibly
    something else. Needless to say, it is a good idea to cross-check
    the definitions when working with theta functions.

    **Definition**

    The four Jacobi theta functions as implemented by :func:`jtheta`
    are defined by the following infinite series:

    .. math ::

      \vartheta_1(z,q) = 2 q^{1/4} \sum_{n=0}^{\infty}
        (-1)^n q^{n^2+n\,} \sin((2n+1)z)

      \vartheta_2(z,q) = 2 q^{1/4} \sum_{n=0}^{\infty}
        q^{n^{2\,} + n} \cos((2n+1)z)

      \vartheta_3(z,q) = 1 + 2 \sum_{n=0}^{\infty}
        q^{n^2\,} \cos(2 n z)

      \vartheta_4(z,q) = 1 + 2 \sum_{n=0}^{\infty}
        (-q)^{n^2\,} \cos(2 n z)

    For `|q| \ll 1`, these series converge very quickly, so the
    Jacobi theta functions can efficiently be evaluated to high
    precision.

    **Examples and basic properties**

    Considered as functions of `z`, the Jacobi theta functions may be
    viewed as generalizations of the ordinary trigonometric functions
    cos and sin. They are periodic functions::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print jtheta(1, 0.1, 1/5.)
        0.117756191842059
        >>> print jtheta(1, 0.1 + 2*pi, 1/5.)
        0.117756191842059

    Indeed, the series defining the theta functions are essentially
    trigonometric Fourier series. The coefficients can be retrieved
    using :func:`fourier`::

        >>> nprint(fourier(lambda x: jtheta(2, x, 0.5), [-pi, pi], 4))
        ([0.0, 1.68179, 0.0, 0.420448, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0])

    The Jacobi theta functions are also so-called quasiperiodic
    functions of `z` and `\tau`, meaning that for fixed `\tau`,
    `\vartheta_n(z, q)` and `\vartheta_n(z+\pi \tau, q)` are the same
    except for an exponential factor::

        >>> tau = 0.3*j
        >>> q = exp(pi*j*tau)
        >>> z = 10
        >>> print jtheta(4, z+tau*pi, q)
        (-0.682420280786035 + 1.5266839997214j)
        >>> print -exp(-2*j*z)/q * jtheta(4, z, q)
        (-0.682420280786035 + 1.5266839997214j)

    The Jacobi theta functions satisfy a huge number of other
    functional equations, such as the following identity (valid for
    any `q`)::

        >>> q = 0.3
        >>> print jtheta(3,0,q)**4
        6.82374408935276
        >>> print jtheta(2,0,q)**4 + jtheta(4,0,q)**4
        6.82374408935276

    Extensive listings of identities satisfied by the Jacobi theta
    functions can be found in standard reference works.

    The Jacobi theta functions are related to the gamma function
    for special arguments::

        >>> print jtheta(3, 0, exp(-pi))
        1.08643481121331
        >>> print pi**(1/4.) / gamma(3/4.)
        1.08643481121331

    :func:`jtheta` supports arbitrary precision evaluation and complex
    arguments::

        >>> mp.dps = 50
        >>> print jtheta(4, sqrt(2), 0.5)
        2.0549510717571539127004115835148878097035750653737
        >>> mp.dps = 25
        >>> print jtheta(4, 1+2j, (1+j)/5)
        (7.180331760146805926356634 - 1.634292858119162417301683j)

    **Possible issues**

    For `|q| \ge 1` or `\Im(\tau) \le 0`, :func:`jtheta` raises
    ``ValueError``. This exception is also raised for `|q|` extremely
    close to 1 (or equivalently `\tau` very close to 0), since the
    series would converge too slowly::

        >>> jtheta(1, 10, 0.99999999 * exp(0.5*j))
        Traceback (most recent call last):
          ...
        ValueError: abs(q) > Q_LIM = 1.000000

    """
    z = mpmathify(z)
    q = mpmathify(q)

    # Implementation note
    # If z.imag is close to zero, _jacobi_theta2 and _jacobi_theta3
    # are used,
    # which compute the series starting from n=0 using fixed precision
    # numbers;
    # otherwise  _jacobi_theta2a and _jacobi_theta3a are used, which compute
    # the series starting from n=n0, which is the largest term.

    # TODO write _jacobi_theta2a and _jacobi_theta3a using fixed precision.

    if abs(q) > Q_LIM:
        raise ValueError('abs(q) > Q_LIM = %f' % Q_LIM)

    extra = 10
    cz = 0.5
    extra2 = 50
    prec0 = mp.prec
    try:
        mp.prec += extra
        if n == 1:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(log(q).real):
                    mp.dps += extra2
                    res = _jacobi_theta2(z - pi/2, q)
                else:
                    mp.dps += 10
                    res = _jacobi_theta2a(z - pi/2, q)
            else:
                res = _jacobi_theta2(z - pi/2, q)
        elif n == 2:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(log(q).real):
                    mp.dps += extra2
                    res = _jacobi_theta2(z, q)
                else:
                    mp.dps += 10
                    res = _jacobi_theta2a(z, q)
            else:
                res = _jacobi_theta2(z, q)
        elif n == 3:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(log(q).real):
                    mp.dps += extra2
                    res = _jacobi_theta3(z, q)
                else:
                    mp.dps += 10
                    res = _jacobi_theta3a(z, q)
            else:
                res = _jacobi_theta3(z, q)
        elif n == 4:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(log(q).real):
                    mp.dps += extra2
                    res = _jacobi_theta3(z, -q)
                else:
                    mp.dps += 10
                    res = _jacobi_theta3a(z, -q)
            else:
                res = _jacobi_theta3(z, -q)
        else:
            raise ValueError
    finally:
        mp.prec = prec0
    return res

def djtheta(n, z, q, nd=1):
    r"""
    For an integer `nd \ge 1`, computes the `nd`:th derivative with
    respect to `z` of the Jacobi theta function `\vartheta_n(z,q)`::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print djtheta(3, 7, 0.2)
        -0.795947847483158
        >>> print diff(lambda x: jtheta(3, x, 0.2), 7)
        -0.795947847483158

    For additional details, see :func:`jtheta`.
    """

    z = mpmathify(z)
    q = mpmathify(q)

    if abs(q) > Q_LIM:
        raise ValueError('abs(q) > Q_LIM = %f' % Q_LIM)
    extra = 10 + mp.prec * nd // 10
    cz = 0.5
    extra2 = 50
    prec0 = mp.prec
    try:
        mp.prec += extra
        if n == 1:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(log(q).real):
                    mp.dps += extra2
                    res = _djacobi_theta2(z - pi/2, q, nd)
                else:
                    mp.dps += 10
                    res = _djacobi_theta2a(z - pi/2, q, nd)
            else:
                res = _djacobi_theta2(z - pi/2, q, nd)
        elif n == 2:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(log(q).real):
                    mp.dps += extra2
                    res = _djacobi_theta2(z, q, nd)
                else:
                    mp.dps += 10
                    res = _djacobi_theta2a(z, q, nd)
            else:
                res = _djacobi_theta2(z, q, nd)
        elif n == 3:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(log(q).real):
                    mp.dps += extra2
                    res = _djacobi_theta3(z, q, nd)
                else:
                    mp.dps += 10
                    res = _djacobi_theta3a(z, q, nd)
            else:
                res = _djacobi_theta3(z, q, nd)
        elif n == 4:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(log(q).real):
                    mp.dps += extra2
                    res = _djacobi_theta3(z, -q, nd)
                else:
                    mp.dps += 10
                    res = _djacobi_theta3a(z, -q, nd)
            else:
                res = _djacobi_theta3(z, -q, nd)
        else:
            raise ValueError
    finally:
        mp.prec = prec0
    return res

def jsn(u, m):
    """
    Computes of the Jacobi elliptic sn function in terms
    of Jacobi theta functions.
    `u` is any complex number, `m` must be in the unit disk

    The sn-function is doubly periodic in the complex
    plane with periods `4 K(m)` and `2 i K(1-m)`
    (see :func:`ellipk`)::

        >>> from mpmath import *
        >>> mp.dps = 25
        >>> print jsn(2, 0.25)
        0.9628981775982774425751399
        >>> print jsn(2+4*ellipk(0.25), 0.25)
        0.9628981775982774425751399
        >>> print chop(jsn(2+2*j*ellipk(1-0.25), 0.25))
        0.9628981775982774425751399

    """
    if abs(m) < eps:
        return sin(u)
    elif m == one:
        return tanh(u)
    else:
        extra = 10
    try:
        mp.prec += extra
        q = calculate_nome(sqrt(m))

        v3 = jtheta(3, zero, q)
        v2 = jtheta(2, zero, q)        # mathworld says v4
        arg1 = u / (v3*v3)
        v1 = jtheta(1, arg1, q)
        v4 = jtheta(4, arg1, q)

        sn = (v3/v2)*(v1/v4)
    finally:
        mp.prec -= extra

    return sn


def jcn(u, m):
    """
    Computes of the Jacobi elliptic cn function in terms
    of Jacobi theta functions.
    `u` is any complex number, `m` must be in the unit disk

    The cn-function is doubly periodic in the complex
    plane with periods `4 K(m)` and `4 i K(1-m)`
    (see :func:`ellipk`)::

        >>> from mpmath import *
        >>> mp.dps = 25
        >>> print jcn(2, 0.25)
        -0.2698649654510865792581416
        >>> print jcn(2+4*ellipk(0.25), 0.25)
        -0.2698649654510865792581416
        >>> print chop(jcn(2+4*j*ellipk(1-0.25), 0.25))
        -0.2698649654510865792581416

    """
    if abs(m) < eps:
        return cos(u)
    elif m == one:
        return sech(u)
    else:
        extra = 10
    try:
        mp.prec += extra
        q = calculate_nome(sqrt(m))

        v3 = jtheta(3, zero, q)
        v2 = jtheta(2, zero, q)
        v04 = jtheta(4, zero, q)

        arg1 = u / (v3*v3)

        v1 = jtheta(2, arg1, q)
        v4 = jtheta(4, arg1, q)

        cn = (v04/v2)*(v1/v4)
    finally:
        mp.prec -= extra
    return cn


def jdn(u, m):
    """
    Computes of the Jacobi elliptic dn function in terms
    of Jacobi theta functions.
    `u` is any complex number, `m` must be in the unit disk

    The dn-function is doubly periodic in the complex
    plane with periods `2 K(m)` and `4 i K(1-m)`
    (see :func:`ellipk`)::

        >>> from mpmath import *
        >>> mp.dps = 25
        >>> print jdn(2, 0.25)
        0.8764740583123262286931578
        >>> print jdn(2+2*ellipk(0.25), 0.25)
        0.8764740583123262286931578
        >>> print chop(jdn(2+4*j*ellipk(1-0.25), 0.25))
        0.8764740583123262286931578

    """
    if m == zero:
        return one
    elif m == one:
        return sech(u)
    else:
        extra = 10
    try:
        mp.prec += extra
        q = calculate_nome(sqrt(m))

        v3 = jtheta(3, zero, q)
        v2 = jtheta(2, zero, q)
        v04 = jtheta(4, zero, q)

        arg1 = u / (v3*v3)

        v1 = jtheta(3, arg1, q)
        v4 = jtheta(4, arg1, q)

        cn = (v04/v3)*(v1/v4)
    finally:
        mp.prec -= extra
    return cn

jacobi_elliptic_sn = jsn
jacobi_elliptic_cn = jcn
jacobi_elliptic_dn = jdn

if __name__ == '__main__':
    import doctest
    doctest.testmod()
