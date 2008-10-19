"""
    elliptic.py

    Implements the Jacobi theta and Jacobi elliptic functions, using
    arbitrary precision math library

    Author: M.T. Taschuk

    References:

    [1] Abramowitz & Stegun. 'Handbook of Mathematical Functions, 9th Ed.',
        (Dover duplicate of 1972 edition)
    [2] Whittaker 'A Course of Modern Analysis, 4th Ed.', 1946,
        Cambridge Univeristy Press

"""
import sys

from mptypes import (mpf, convert_lossless, eps)
from functions import (pi, sqrt, cos, sin, exp, ellipk, sech)

def calculate_nome(k):
    """
    Calculate the nome, q, from the value for k.

    Useful factoids:

    k**2 = m;   m is used in Abramowitz
    """
    k = convert_lossless(k)

    if k > mpf('1'):             # range error
        raise ValueError

    zero = mpf('0')
    one = mpf('1')

    if k == zero:
        return zero
    elif k == one:
        return one
    else:
        kprimesquared = one - k**2
        kprime = sqrt(kprimesquared)
        top = ellipk(kprimesquared)
        bottom = ellipk(k**2)

        argument = mpf('-1')*pi*top/bottom

        nome = exp(argument)
        return nome

def calculate_k(q, verbose=False):
    """
    Calculates the value of k for a particular nome, q.

    Uses special cases of the jacobi theta functions, with
    q as an argument, rather than m.

    k = (v2(0, q)/v3(0, q))**2
    """
    zero = mpf('0')
    one = mpf('1')

    q = convert_lossless(q)
    if q > one or q < zero:
        raise ValueError

    # calculate v2(0, q)
    sum = zero
    term = zero                     # series starts at zero
    while True:
        factor1 = q**(term*(term + 1))
        term_n = factor1            # suboptimal, kept for readability
        sum = sum + term_n

        if verbose:
            print >> sys.stderr, '\tTerm: %d' % term,
            print >> sys.stderr, '\tterm_n: %e' % term_n,
            print >> sys.stderr, '\tsum: %e' % sum

        if factor1 == zero:         # all further terms will be zero
            break
        #if log(term_n, '10') < -1*mpf.dps:
        if abs(term_n) < eps:
            break

        term = term + 1

    v2 = 2*q**(mpf('0.25'))*sum

    # calculate v3(0, q)
    sum = zero
    term = one                          # series starts at one
    while True:
        factor1 = q**(term*term)
        term_n = factor1                # suboptimal, kept for readability
        sum = sum + term_n

        if factor1 == mpf('0'):      # all further terms will be zero
            break
        #if log(term_n, '10') < -1*mpf.dps:
        if abs(term_n) < eps:
            break

        term = term + 1

    v3 = one + 2*sum

    k = v2**2/v3**2

    return k

def jacobi_theta_1(z, m, verbose=False):
    """
    Implements the jacobi theta function 1, using the series expansion
    found in Abramowitz & Stegun [1].

    z is any complex number, but only reals here?
    m is the parameter, which must be converted to the nome
    """
    if verbose:
        print >> sys.stderr, 'elliptic.jacobi_theta_1'

    m = convert_lossless(m)
    z = convert_lossless(z)

    k = sqrt(m)
    q = calculate_nome(k)

    if verbose:
        print >> sys.stderr, '\tk: %f ' % k
        print >> sys.stderr, '\tq: %f ' % q

    if abs(q) >= mpf('1'):
        raise ValueError

    sum = 0
    term = 0                                    # series starts at 0

    minusone = mpf('-1')
    zero = mpf('0')
    one = mpf('1')

    if z == zero:
        if verbose:
            print >> sys.stderr, 'elliptic.jacobi_theta_1: z == 0, return 0'
        return zero
    elif q == zero:
        if verbose:
            print >> sys.stderr, 'elliptic.jacobi_theta_1: q == 0, return 0'
        return zero
    else:
        if verbose:
            print >> sys.stderr, 'ellipticl.jacobi_theta_1: calculating'
        while True:
            if (term % 2) == 0:
                factor0 = one
            else:
                factor0 = minusone

            factor1 = q**(term*(term +1))
            factor2 = sin((2*term + 1)*z)

            term_n = factor0 * factor1 * factor2
            sum = sum + term_n
            if verbose:
                print >> sys.stderr, '\tTerm: %d' % term,
                print >> sys.stderr, '\tterm_n: %e' % term_n,
                print >> sys.stderr, '\tsum: %e' % sum

            if factor1 == zero:         # all further terms will be zero
                break
            if factor2 != zero:         # check precision iff cos != 0
                #if log(term_n, '10') < -1*mpf.dps:
                if abs(term_n) < eps:
                    break

            term = term + 1

        return (2*q**(0.25))*sum

    # can't get here
    print >> sys.stderr, 'elliptic.jacobi_theta_1 in impossible state'
    sys.exit(2)

def jacobi_theta_2(z, m, verbose=False):
    """
    Implements the jacobi theta function 2, using the series expansion
    found in Abramowitz & Stegun [4].

    z is any complex number, but only reals here?
    m is the parameter, which must be converted to the nome
    """
    if verbose:
        print >> sys.stderr, 'elliptic.jacobi_theta_2'

    m = convert_lossless(m)
    z = convert_lossless(z)

    k = sqrt(m)
    q = calculate_nome(k)

    if verbose:
        print >> sys.stderr, '\tk: %f ' % k
        print >> sys.stderr, '\tq: %f ' % q

    if abs(q) >= mpf('1'):
        raise ValueError

    sum = 0
    term = 0                    # series starts at zero

    zero = mpf('0')

    if q == zero:
        if verbose:
            print >> sys.stderr, 'elliptic.jacobi_theta_2: q == 0, return 0'
        return zero
    else:
        if verbose:
            print >> sys.stderr, 'elliptic.jacobi_theta_2: calculating'
        if z == zero:
            factor2 = mpf('1')
            while True:
                factor1 = q**(term*(term + 1))
                term_n = factor1            # suboptimal, kept for readability
                sum = sum + term_n

                if verbose:
                    print >> sys.stderr, '\tTerm: %d' % term,
                    print >> sys.stderr, '\tterm_n: %e' % term_n,
                    print >> sys.stderr, '\tsum: %e' % sum

                if factor1 == zero:         # all further terms will be zero
                    break
                #if log(term_n, '10') < -1*mpf.dps:
                if abs(term_n) < eps:
                    break

                term = term + 1
        else:
            while True:
                factor1 = q**(term*(term + 1))
                if z == zero:
                    factor2 = 1
                else:
                    factor2 = cos((2*term + 1)*z)

                term_n = factor1 * factor2
                sum = sum + term_n

                if verbose:
                    print >> sys.stderr, '\tTerm: %d' % term,
                    print >> sys.stderr, '\tterm_n: %e' % term_n,
                    print >> sys.stderr, '\tsum: %e' % sum

                if factor1 == zero:         # all further terms will be zero
                    break
                if factor2 != zero:         # check precision iff cos != 0
                    #if log(term_n, '10') < -1*mpf.dps:
                    if abs(term_n) < eps:
                        break

                term = term + 1

        return (2*q**(0.25))*sum

    # can't get here
    print >> sys.stderr, 'elliptic.jacobi_theta_2 in impossible state'
    sys.exit(2)

def jacobi_theta_3(z, m):
    """
    Implements the jacobi theta function 2, using the series expansion
    found in Abramowitz & Stegun [4].

    z is any complex number, but only reals here?
    m is the parameter, which must be converted to the nome
    """
    m = convert_lossless(m)
    z = convert_lossless(z)

    k = sqrt(m)
    q = calculate_nome(k)

    if abs(q) >= mpf('1'):
        raise ValueError

    delta_sum = 1
    term = 1                    # series starts at 1

    zero = mpf('0')
    sum = zero

    if z == zero:
        factor2 = mpf('1')
        while True:

            factor1 = q**(term*term)
            term_n = factor1                # suboptimal, kept for readability
            sum = sum + term_n

            if factor1 == mpf('0'):      # all further terms will be zero
                break
            #if log(term_n, '10') < -1*mpf.dps:
            if abs(term_n) < eps:
                break

            term = term + 1

    else:
        while True:
            factor1 = q**(term*term)
            if z == zero:
                factor2 = 1
            else:
                factor2 = cos(2*term*z)

            term_n = factor1 * factor2
            sum = sum + term_n

            if factor1 == mpf('0'):      # all further terms will be zero
                break
            if factor2 != mpf('0'):      # check precision iff cos != 0
                #if log(term_n, '10') < -1*mpf.dps:
                if abs(term_n) < eps:
                    break

            term = term + 1

    return 1 + 2*sum

def jacobi_theta_4(z, m):
    """
    Implements the series expansion of the jacobi theta function
    1, where z == 0.

    z is any complex number, but only reals here?
    m is the parameter, which must be converted to the nome
    """
    m = convert_lossless(m)
    z = convert_lossless(z)

    k = sqrt(m)
    q = calculate_nome(k)

    if abs(q) >= mpf('1'):
        raise ValueError

    sum = 0
    delta_sum = 1
    term = 1                    # series starts at 1

    zero = mpf('0')

    if z == zero:
        factor2 = mpf('1')
        while True:
            if (term % 2) == 0:
                factor0 = 1
            else:
                factor0 = -1
            factor1 = q**(term*term)
            term_n = factor0 * factor1 * factor2
            sum = sum + term_n

            if factor1 == mpf('0'):      # all further terms will be zero
                break
            #if log(term_n, '10') < -1*mpf.dps:
            if abs(term_n) < eps:
                break

            term = term + 1
    else:
        while True:
            if (term % 2) == 0:
                factor0 = 1
            else:
                factor0 = -1
            factor1 = q**(term*term)
            if z == zero:
                factor2 = 1
            else:
                factor2 = cos(2*term*z)

            term_n = factor0 * factor1 * factor2
            sum = sum + term_n

            if factor1 == mpf('0'):      # all further terms will be zero
                break
            if factor2 != mpf('0'):      # check precision iff cos != 0
                #if log(term_n, '10') < -1*mpf.dps:
                if abs(term_n) < eps:
                    break

            term = term + 1

    return 1 + 2*sum

def jacobi_elliptic_sn(u, m, verbose=False):
    """
    Implements the jacobi elliptic sn function, using the expansion in
    terms of q, from Abramowitz 16.23.1.

    u is any complex number, m is the parameter

    Alternative implementation:

    Expansion in terms of jacobi theta functions appears to fail with
    round off error, despite   I also think that the expansion in
    terms of q is much faster than four expansions in terms of q.

    **********************************
    Previous implementation kept here:

    if not isinstance(u, mpf):
        raise TypeError
    if not isinstance(m, mpf):
        raise TypeError

    zero = mpf('0')

    if u == zero and m == 0:
        return zero
    else:
        q = calculate_nome(sqrt(m))

        v3 = jacobi_theta_3(zero, q)
        v2 = jacobi_theta_2(zero, q)        # mathworld says v4

        arg1 = u / (v3*v3)

        v1 = jacobi_theta_1(arg1, q)
        v4 = jacobi_theta_4(arg1, q)

        sn = (v3/v2)*(v1/v4)

        return sn
    **********************************
    """
    u = convert_lossless(u)
    m = convert_lossless(m)

    if verbose:
        print >> sys.stderr, '\nelliptic.jacobi_elliptic_sn'
        print >> sys.stderr, '\tu: %1.12f' % u
        print >> sys.stderr, '\tm: %1.12f' % m

    zero = mpf('0')
    onehalf = mpf('0.5')
    one = mpf('1')
    two = mpf('2')

    if m == zero:           # sn collapes to sin(u)
        if verbose:
            print >> sys.stderr, '\nsn: special case, m == 0'
        return sin(u)
    elif m == one:          # sn collapses to tanh(u)
        if verbose:
            print >> sys.stderr, '\nsn: special case, m == 1'
        return tanh(u)
    else:
        k = sqrt(m)                        # convert m to k
        q = calculate_nome(k)
        v = (pi * u) / (two*ellipk(k**2))

    if v == pi or v == zero:     # sin factor always zero
        return zero

    sum = zero
    term = zero                     # series starts at zero

    while True:
        if verbose:
            print >> sys.stderr, 'elliptic.jacobi_elliptic_sn: calculating'
        factor1 = (q**(term + onehalf)) / (one - q**(two*term + one))
        factor2 = sin((two*term + one)*v)

        term_n = factor1*factor2
        sum = sum + term_n

        if verbose:
            print >> sys.stderr, '\tTerm: %d' % term,
            print >> sys.stderr, '\tterm_n: %e' % term_n,
            print >> sys.stderr, '\tsum: %e' % sum

        if not factor2 == zero:
            #if log(term_n, '10') < -1*mpf.dps:
            if abs(term_n) < eps:
                break

        term = term + one

    answer = (two*pi) / (sqrt(m) * ellipk(k**2)) * sum

    return answer

def jacobi_elliptic_cn(u, m, verbose=False):
    """
    Implements the jacobi elliptic cn function, using the expansion in
    terms of q, from Abramowitz 16.23.2.
    """
    u = convert_lossless(u)
    m = convert_lossless(m)

    if verbose:
        print >> sys.stderr, '\nelliptic.jacobi_elliptic_cn'
        print >> sys.stderr, '\tu: %1.12f' % u
        print >> sys.stderr, '\tm: %1.12f' % m

    zero = mpf('0')
    onehalf = mpf('0.5')
    one = mpf('1')
    two = mpf('2')

    if m == zero:                   # cn collapses to cos(u)
        if verbose:
            print >> sys.stderr, 'cn: special case, m == 0'
        return cos(u)
    elif m == one:                  # cn collapses to sech(u)
        if verbose:
            print >> sys.stderr, 'cn: special case, m == 1'
        return sech(u)
    else:
        k = sqrt(m)                        # convert m to k
        q = calculate_nome(k)
        kprimesquared = one - k**2
        kprime = sqrt(kprimesquared)
        v = (pi * u) / (two*ellipk(k**2))

    sum = zero
    term = zero                     # series starts at zero

    if verbose:
        print >> sys.stderr, 'elliptic.jacobi_elliptic_cn: calculating'
    while True:
        factor1 = (q**(term + onehalf)) / (one + q**(two*term + one))
        factor2 = cos((two*term + one)*v)

        term_n = factor1*factor2
        sum = sum + term_n

        if verbose:
            print >> sys.stderr, '\tTerm: %d' % term,
            print >> sys.stderr, '\tterm_n: %e' % term_n,
            print >> sys.stderr, '\tsum: %e' % sum

        if not factor2 == zero:
            #if log(term_n, '10') < -1*mpf.dps:
            if abs(term_n) < eps:
                break

        term = term + one

    answer = (two*pi) / (sqrt(m) * ellipk(k**2)) * sum

    return answer

def jacobi_elliptic_dn(u, m, verbose=False):
    """
    Implements the jacobi elliptic cn function, using the expansion in
    terms of q, from Abramowitz 16.23.3.
    """
    u = convert_lossless(u)
    m = convert_lossless(m)

    if verbose:
        print >> sys.stderr, '\nelliptic.jacobi_elliptic_dn'
        print >> sys.stderr, '\tu: %1.12f' % u
        print >> sys.stderr, '\tm: %1.12f' % m

    zero = mpf('0')
    onehalf = mpf('0.5')
    one = mpf('1')
    two = mpf('2')

    if m == zero:           # dn collapes to 1
        return one
    elif m == one:          # dn collapses to sech(u)
        return sech(u)
    else:
        k = sqrt(m)                        # convert m to k
        q = calculate_nome(k)
        v = (pi * u) / (two*ellipk(k**2))

    sum = zero
    term = one                  # series starts at one

    if verbose:
        print >> sys.stderr, 'elliptic.jacobi_elliptic_dn: calculating'
    while True:
        factor1 = (q**term) / (one + q**(two*term))
        factor2 = cos(two*term*v)

        term_n = factor1*factor2
        sum = sum + term_n

        if verbose:
            print >> sys.stderr, '\tTerm: %d' % term,
            print >> sys.stderr, '\tterm_n: %e' % term_n,
            print >> sys.stderr, '\tsum: %e' % sum

        if not factor2 == zero:
            #if log(term_n, '10') < -1*mpf.dps:
            if abs(term_n) < eps:
                break

        term = term + one

    K = ellipk(k**2)
    answer = (pi / (two*K)) + (two*pi*sum)/(ellipk(k**2))

    return answer
