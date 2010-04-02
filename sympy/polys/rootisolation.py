"""Tools for real and complex root isolation and refinement. """

from sympy.polys.densebasic import (
    dup_LC, dup_TC, dup_degree,
    dup_strip, dup_reverse,
    dup_convert, dmp_convert,
    dup_terms_gcd,
)

from sympy.polys.densearith import (
    dup_neg, dup_rshift,
)

from sympy.polys.densetools import (
    dup_scale, dup_taylor, dup_mirror,
    dup_transform,
    dup_eval, dmp_eval_in,
    dup_sign_variations,
    dup_clear_denoms,
    dup_real_imag,
    dup_sqf_list,
)

from sympy.polys.factortools import (
    dup_factor_list,
)

from sympy.polys.polyerrors import (
    RefinementFailed,
    DomainError,
)

import operator

def dup_root_upper_bound(f, K):
    """Compute LMQ upper bound for `f`'s positive roots. """
    n, t, P = len(f), K.one, []

    if dup_LC(f, K) < 0:
        f = dup_neg(f, K)

    f = list(reversed(f))

    for i in xrange(0, n):
        if f[i] >= 0:
            continue

        a, Q = K.log(-f[i], 2), []

        for j in xrange(i+1, n):

            if f[j] <= 0:
                continue

            q = t + a - K.log(f[j], 2)
            Q.append(q // (j - i))

            t += 1

        if not Q:
            continue

        P.append(min(Q))

    if not P:
        return None
    else:
        return 2.0**(max(P)+1)

def dup_root_lower_bound(f, K):
    """Compute LMQ lower bound for `f`'s positive roots. """
    bound = dup_root_upper_bound(dup_reverse(f), K)

    if bound is not None:
        return 1.0 / bound
    else:
        return None

def _mobius_from_interval(I, field):
    """Convert an open interval to a Mobius transform. """
    s, t = I

    a, c = field.numer(s), field.denom(s)
    b, d = field.numer(t), field.denom(t)

    return a, b, c, d

def _mobius_to_interval(M, field):
    """Convert a Mobius transform to an open interval. """
    a, b, c, d = M

    s, t = field(a, c), field(b, d)

    if s <= t:
        return (s, t)
    else:
        return (t, s)

def dup_step_refine_real_root(f, M, K, fast=False):
    """One step of positive real root refinement algorithm. """
    a, b, c, d = M

    if a == b and c == d:
        return f, (a, b, c, d)

    A = dup_root_lower_bound(f, K)

    if A is not None:
        A = K(int(A))
    else:
        A = K.zero

    if fast and A > 16:
        f = dup_scale(f, A, K)
        a, c, A = A*a, A*c, K.one

    if A >= K.one:
        f = dup_taylor(f, A, K)
        b, d = A*a + b, A*c + d

        if not dup_eval(f, K.zero, K):
            return f, (b, b, d, d)

    f, g = dup_taylor(f, K.one, K), f

    a1, b1, c1, d1 = a, a+b, c, c+d

    if not dup_eval(f, K.zero, K):
        return f, (b1, b1, d1, d1)

    k = dup_sign_variations(f, K)

    if k == 1:
        a, b, c, d = a1, b1, c1, d1
    else:
        f = dup_taylor(dup_reverse(g), K.one, K)

        if not dup_eval(f, K.zero, K):
            f = dup_rshift(f, 1, K)

        a, b, c, d = b, a+b, d, c+d

    return f, (a, b, c, d)

def dup_inner_refine_real_root(f, M, K, eps=None, steps=None, disjoint=None, fast=False, mobius=False):
    """Refine a positive root of `f` given a Mobius transform or an interval. """
    F = K.get_field()

    if len(M) == 2:
        a, b, c, d = _mobius_from_interval(M, F)
    else:
        a, b, c, d = M

    while not c:
        f, (a, b, c, d) = dup_step_refine_real_root(f, (a, b, c, d), K, fast=fast)

    if eps is not None and steps is not None:
        for i in xrange(0, steps):
            if abs(F(a, c) - F(b, d)) >= eps:
                f, (a, b, c, d) = dup_step_refine_real_root(f, (a, b, c, d), K, fast=fast)
            else:
                break
    else:
        if eps is not None:
            while abs(F(a, c) - F(b, d)) >= eps:
                f, (a, b, c, d) = dup_step_refine_real_root(f, (a, b, c, d), K, fast=fast)

        if steps is not None:
            for i in xrange(0, steps):
                f, (a, b, c, d) = dup_step_refine_real_root(f, (a, b, c, d), K, fast=fast)

    if disjoint is not None:
        while True:
            u, v = _mobius_to_interval((a, b, c, d), F)

            if v <= disjoint or disjoint <= u:
                break
            else:
                f, (a, b, c, d) = dup_step_refine_real_root(f, (a, b, c, d), K, fast=fast)

    if not mobius:
        return _mobius_to_interval((a, b, c, d), F)
    else:
        return f, (a, b, c, d)

def dup_outer_refine_real_root(f, s, t, K, eps=None, steps=None, disjoint=None, fast=False):
    """Refine a positive root of `f` given an interval `(s, t)`. """
    a, b, c, d = _mobius_from_interval((s, t), K.get_field())

    f = dup_transform(f, dup_strip([a, b]),
                         dup_strip([c, d]), K)

    if dup_sign_variations(f, K) != 1:
        raise RefinementFailed("there should be exactly one root in (%s, %s) interval" % (s, t))

    return dup_inner_refine_real_root(f, (a, b, c, d), K, eps=eps, steps=steps, disjoint=disjoint, fast=fast)

def dup_refine_real_root(f, s, t, K, eps=None, steps=None, disjoint=None, fast=False):
    """Refine real root's approximating interval to the given precision. """
    if K.is_QQ:
        (_, f), K = dup_clear_denoms(f, K, convert=True), K.get_ring()
    elif not K.is_ZZ:
        raise DomainError("real root refinement not supported over %s" % K)

    if s == t:
        return (s, t)

    if s > t:
        s, t = t, s

    negative = False

    if s < 0:
        if t <= 0:
            f, s, t, negative = dup_mirror(f, K), -t, -s, True
        else:
            raise ValueError("can't refine a real root in (%s, %s)" % (s, t))

    if negative and disjoint is not None:
        if disjoint < 0:
            disjoint = -disjoint
        else:
            disjoint = None

    s, t = dup_outer_refine_real_root(f, s, t, K, eps=eps, steps=steps, disjoint=disjoint, fast=fast)

    if negative:
        return (-t, -s)
    else:
        return ( s,  t)

def dup_inner_isolate_real_roots(f, K, eps=None, fast=False):
    """Internal function for isolation positive roots up to given precision. """
    a, b, c, d = K.one, K.zero, K.zero, K.one

    k = dup_sign_variations(f, K)

    if k == 0:
        return []
    if k == 1:
        roots = [dup_inner_refine_real_root(f, (a, b, c, d), K, eps=eps, fast=fast, mobius=True)]
    else:
        roots, stack = [], [(a, b, c, d, f, k)]

        F = K.get_field()

        while stack:
            a, b, c, d, f, k = stack.pop()

            A = dup_root_lower_bound(f, K)

            if A is not None:
                A = K(int(A))
            else:
                A = K.zero

            if fast and A > 16:
                f = dup_scale(f, A, K)
                a, c, A = A*a, A*c, K.one

            if A >= K.one:
                f = dup_taylor(f, A, K)
                b, d = A*a + b, A*c + d

                if not dup_TC(f, K):
                    roots.append((f, (b, b, d, d)))
                    f = dup_rshift(f, 1, K)

                k = dup_sign_variations(f, K)

                if k == 0:
                    continue
                if k == 1:
                    roots.append(dup_inner_refine_real_root(f, (a, b, c, d), K, eps=eps, fast=fast, mobius=True))
                    continue

            f1 = dup_taylor(f, K.one, K)

            a1, b1, c1, d1, r = a, a+b, c, c+d, 0

            if not dup_TC(f1, K):
                roots.append((f1, (b1, b1, d1, d1)))
                f1, r = dup_rshift(f1, 1, K), 1

            k1 = dup_sign_variations(f1, K)
            k2 = k - k1 - r

            a2, b2, c2, d2 = b, a+b, d, c+d

            if k2 > 1:
                f2 = dup_taylor(dup_reverse(f), K.one, K)

                if not dup_TC(f2, K):
                    f2 = dup_rshift(f2, 1, K)

                k2 = dup_sign_variations(f2, K)
            else:
                f2 = None

            if k1 < k2:
                a1, a2, b1, b2 = a2, a1, b2, b1
                c1, c2, d1, d2 = c2, c1, d2, d1
                f1, f2, k1, k2 = f2, f1, k2, k1

            if not k1:
                continue

            if f1 is None:
                f1 = dup_taylor(dup_reverse(f), K.one, K)

                if not dup_TC(f1, K):
                    f1 = dup_rshift(f1, 1, K)

            if k1 == 1:
                roots.append(dup_inner_refine_real_root(f1, (a1, b1, c1, d1), K, eps=eps, fast=fast, mobius=True))
            else:
                stack.append((a1, b1, c1, d1, f1, k1))

            if not k2:
                continue

            if f2 is None:
                f2 = dup_taylor(dup_reverse(f), K.one, K)

                if not dup_TC(f2, K):
                    f2 = dup_rshift(f2, 1, K)

            if k2 == 1:
                roots.append(dup_inner_refine_real_root(f2, (a2, b2, c2, d2), K, eps=eps, fast=fast, mobius=True))
            else:
                stack.append((a2, b2, c2, d2, f2, k2))

    return roots

def _discard_if_outside_interval(f, M, inf, sup, K, negative, fast, mobius):
    """Discard an isolating interval if outside ``(inf, sup)``. """
    F = K.get_field()

    while True:
        u, v = _mobius_to_interval(M, F)

        if negative:
            u, v = -v, -u

        if (inf is None or u >= inf) and (sup is None or v <= sup):
            if not mobius:
                return u, v
            else:
                return f, M
        elif (sup is not None and u > sup) or (inf is not None and v < inf):
            return None
        else:
            f, M = dup_step_refine_real_root(f, M, K, fast=fast)

def dup_inner_isolate_positive_roots(f, K, eps=None, inf=None, sup=None, fast=False, mobius=False):
    """Iteratively compute disjoint positive root isolation intervals. """
    if sup is not None and sup < 0:
        return []

    roots = dup_inner_isolate_real_roots(f, K, eps=eps, fast=fast)

    F, results = K.get_field(), []

    if inf is not None or sup is not None:
        for f, M in roots:
            result = _discard_if_outside_interval(f, M, inf, sup, K, False, fast, mobius)

            if result is not None:
                results.append(result)
    elif not mobius:
        for f, M in roots:
            u, v = _mobius_to_interval(M, F)
            results.append((u, v))
    else:
        results = roots

    return results

def dup_inner_isolate_negative_roots(f, K, inf=None, sup=None, eps=None, fast=False, mobius=False):
    """Iteratively compute disjoint negative root isolation intervals. """
    if inf is not None and inf >= 0:
        return []

    roots = dup_inner_isolate_real_roots(dup_mirror(f, K), K, eps=eps, fast=fast)

    F, results = K.get_field(), []

    if inf is not None or sup is not None:
        for f, M in roots:
            result = _discard_if_outside_interval(f, M, inf, sup, K, True, fast, mobius)

            if result is not None:
                results.append(result)
    elif not mobius:
        for f, M in roots:
            u, v = _mobius_to_interval(M, F)
            results.append((-v, -u))
    else:
        results = roots

    return results

def _isolate_zero(f, K, inf, sup, basis=False, sqf=False):
    """Handle special case of CF algorithm when ``f`` is homogeneous. """
    j, f = dup_terms_gcd(f, K)

    if j > 0:
        F = K.get_field()

        if (inf is None or inf <= 0) and (sup is None or 0 <= sup):
            if not sqf:
                if not basis:
                    return [((F.zero, F.zero), j)], f
                else:
                    return [((F.zero, F.zero), j, [K.one, K.zero])], f
            else:
                return [(F.zero, F.zero)], f

    return [], f

def dup_isolate_real_roots_sqf(f, K, eps=None, inf=None, sup=None, fast=False):
    """Isolate real roots of a square-free polynomial using CF approach. """
    if K.is_QQ:
        (_, f), K = dup_clear_denoms(f, K, convert=True), K.get_ring()
    elif not K.is_ZZ:
        raise DomainError("isolation of real roots not supported over %s" % K)

    if dup_degree(f) <= 0:
        return []

    I_zero, f = _isolate_zero(f, K, inf, sup, basis=False, sqf=True)

    I_neg = dup_inner_isolate_negative_roots(f, K, eps=eps, inf=inf, sup=sup, fast=fast)
    I_pos = dup_inner_isolate_positive_roots(f, K, eps=eps, inf=inf, sup=sup, fast=fast)

    return sorted(I_neg + I_zero + I_pos)

def dup_isolate_real_roots(f, K, eps=None, inf=None, sup=None, basis=False, fast=False):
    """Isolate real roots using continued fractions approach. """
    if K.is_QQ:
        (_, f), K = dup_clear_denoms(f, K, convert=True), K.get_ring()
    elif not K.is_ZZ:
        raise DomainError("isolation of real roots not supported over %s" % K)

    if dup_degree(f) <= 0:
        return []

    I_zero, f = _isolate_zero(f, K, inf, sup, basis=basis, sqf=False)

    _, factors = dup_sqf_list(f, K)

    if len(factors) == 1:
        ((f, k),) = factors

        I_neg = dup_inner_isolate_negative_roots(f, K, eps=eps, inf=inf, sup=sup, fast=fast)
        I_pos = dup_inner_isolate_positive_roots(f, K, eps=eps, inf=inf, sup=sup, fast=fast)

        I_neg = [ ((u, v), k) for u, v in I_neg ]
        I_pos = [ ((u, v), k) for u, v in I_pos ]
    else:
        I_neg, I_pos = _real_isolate_and_disjoin(factors, K,
            eps=eps, inf=inf, sup=sup, basis=basis, fast=fast)

    return sorted(I_neg + I_zero + I_pos)

def dup_isolate_real_roots_list(polys, K, eps=None, inf=None, sup=None, strict=False, basis=False, fast=False):
    """Isolate real roots of a list of square-free polynomial using CF approach. """
    if K.is_QQ:
        K, F, polys = K.get_ring(), K, polys[:]

        for i, p in enumerate(polys):
            polys[i] = dup_clear_denoms(p, F, K, convert=True)[1]
    elif not K.is_ZZ:
        raise DomainError("isolation of real roots not supported over %s" % K)

    zeros, factors_dict = False, {}

    if (inf is None or inf <= 0) and (sup is None or 0 <= sup):
        zeros, zero_indices = True, {}

    for i, p in enumerate(polys):
        j, p = dup_terms_gcd(p, K)

        if zeros and j > 0:
            zero_indices[i] = j

        for f, k in dup_factor_list(p, K)[1]:
            f = tuple(f)

            if f not in factors_dict:
                factors_dict[f] = {i: k}
            else:
                factors_dict[f][i] = k

    factors_list = []

    for f, indices in factors_dict.items():
        factors_list.append((list(f), indices))

    I_neg, I_pos = _real_isolate_and_disjoin(factors_list, K, eps=eps,
        inf=inf, sup=sup, strict=strict, basis=basis, fast=fast)

    F = K.get_field()

    if not zeros or not zero_indices:
        I_zero = []
    else:
        if not basis:
            I_zero = [((F.zero, F.zero), zero_indices)]
        else:
            I_zero = [((F.zero, F.zero), zero_indices, [K.one, K.zero])]

    return sorted(I_neg + I_zero + I_pos)

def _disjoint_p(M, N, strict=False):
    """Check if Mobius transforms define disjoint intervals. """
    a1, b1, c1, d1 = M
    a2, b2, c2, d2 = N

    a1d1, b1c1 = a1*d1, b1*c1
    a2d2, b2c2 = a2*d2, b2*c2

    if a1d1 == b1c1 and a2d2 == b2c2:
        return True

    if a1d1 > b1c1:
        a1, c1, b1, d1 = b1, d1, a1, c1

    if a2d2 > b2c2:
        a2, c2, b2, d2 = b2, d2, a2, c2

    if not strict:
        return a2*d1 >= c2*b1 or b2*c1 <= d2*a1
    else:
        return a2*d1 >  c2*b1 or b2*c1 <  d2*a1

def _real_isolate_and_disjoin(factors, K, eps=None, inf=None, sup=None, strict=False, basis=False, fast=False):
    """Isolate real roots of a list of polynomials and disjoin intervals. """
    I_pos, I_neg = [], []

    for i, (f, k) in enumerate(factors):
        for F, M in dup_inner_isolate_positive_roots(f, K,
                eps=eps, inf=inf, sup=sup, fast=fast, mobius=True):
            I_pos.append((F, M, k, f))

        for G, N in dup_inner_isolate_negative_roots(f, K,
                eps=eps, inf=inf, sup=sup, fast=fast, mobius=True):
            I_neg.append((G, N, k, f))

    for i, (f, M, k, F) in enumerate(I_pos):
        for j, (g, N, m, G) in enumerate(I_pos[i+1:]):
            while not _disjoint_p(M, N, strict=strict):
                f, M = dup_inner_refine_real_root(f, M, K, steps=1, fast=fast, mobius=True)
                g, N = dup_inner_refine_real_root(g, N, K, steps=1, fast=fast, mobius=True)

            I_pos[i+j+1] = (g, N, m, G)

        I_pos[i] = (f, M, k, F)

    for i, (f, M, k, F) in enumerate(I_neg):
        for j, (g, N, m, G) in enumerate(I_neg[i+1:]):
            while not _disjoint_p(M, N, strict=strict):
                f, M = dup_inner_refine_real_root(f, M, K, steps=1, fast=fast, mobius=True)
                g, N = dup_inner_refine_real_root(g, N, K, steps=1, fast=fast, mobius=True)

            I_neg[i+j+1] = (g, N, m, G)

        I_neg[i] = (f, M, k, F)

    if strict:
        for i, (f, M, k, F) in enumerate(I_neg):
            if not M[0]:
                while not M[0]:
                    f, M = dup_inner_refine_real_root(f, M, K, steps=1, fast=fast, mobius=True)

                I_neg[i] = (f, M, k, F)
                break

        for j, (g, N, m, G) in enumerate(I_pos):
            if not N[0]:
                while not N[0]:
                    g, N = dup_inner_refine_real_root(g, N, K, steps=1, fast=fast, mobius=True)

                I_pos[j] = (g, N, m, G)
                break

    field = K.get_field()

    I_neg = [ (_mobius_to_interval(M, field), k, f) for (_, M, k, f) in I_neg ]
    I_pos = [ (_mobius_to_interval(M, field), k, f) for (_, M, k, f) in I_pos ]

    if not basis:
        I_neg = [ ((-v, -u), k) for ((u, v), k, _) in I_neg ]
        I_pos = [ (( u,  v), k) for ((u, v), k, _) in I_pos ]
    else:
        I_neg = [ ((-v, -u), k, f) for ((u, v), k, f) in I_neg ]
        I_pos = [ (( u,  v), k, f) for ((u, v), k, f) in I_pos ]

    return I_neg, I_pos

OO = 'OO' # Origin of (re, im) coordinate system

Q1 = 'Q1' # Quadrant #1 (++): re > 0 and im > 0
Q2 = 'Q2' # Quadrant #2 (-+): re < 0 and im > 0
Q3 = 'Q3' # Quadrant #3 (--): re < 0 and im < 0
Q4 = 'Q4' # Quadrant #4 (+-): re > 0 and im < 0

A1 = 'A1' # Axis #1 (+0): re > 0 and im = 0
A2 = 'A2' # Axis #2 (0+): re = 0 and im > 0
A3 = 'A3' # Axis #3 (-0): re < 0 and im = 0
A4 = 'A4' # Axis #4 (0-): re = 0 and im < 0

_rules_simple = {
    # A -- CCW --> Q => +1/4 (CCW)
    (A1, Q1): 1,
    (A2, Q2): 1,
    (A3, Q3): 1,
    (A4, Q4): 1,

    # A --  CW --> Q => -1/4 (CCW)
    (A1, Q4): 2,
    (A2, Q1): 2,
    (A3, Q2): 2,
    (A4, Q3): 2,

    # Q -- CCW --> A => +1/4 (CCW)
    (Q1, A2): 3,
    (Q2, A3): 3,
    (Q3, A4): 3,
    (Q4, A1): 3,

    # Q --  CW --> A => -1/4 (CCW)
    (Q1, A1): 4,
    (Q2, A2): 4,
    (Q3, A3): 4,
    (Q4, A4): 4,

    # Q -- CCW --> Q => +1/2 (CCW)
    (Q1, Q2): +5,
    (Q2, Q3): +5,
    (Q3, Q4): +5,
    (Q4, Q1): +5,

    # Q --  CW --> Q => -1/2 (CW)
    (Q1, Q4): -5,
    (Q2, Q1): -5,
    (Q3, Q2): -5,
    (Q4, Q3): -5,
}

_rules_ambiguous = {
    # A -- CCW --> Q => { +1/4 (CCW), -9/4 (CW) }
    (A1, OO, Q1): -1,
    (A2, OO, Q2): -1,
    (A3, OO, Q3): -1,
    (A4, OO, Q4): -1,

    # A --  CW --> Q => { -1/4 (CCW), +7/4 (CW) }
    (A1, OO, Q4): -2,
    (A2, OO, Q1): -2,
    (A3, OO, Q2): -2,
    (A4, OO, Q3): -2,

    # Q -- CCW --> A => { +1/4 (CCW), -9/4 (CW) }
    (Q1, OO, A2): -3,
    (Q2, OO, A3): -3,
    (Q3, OO, A4): -3,
    (Q4, OO, A1): -3,

    # Q --  CW --> A => { -1/4 (CCW), +7/4 (CW) }
    (Q1, OO, A1): -4,
    (Q2, OO, A2): -4,
    (Q3, OO, A3): -4,
    (Q4, OO, A4): -4,

    # A --  OO --> A => { +1/1 (CCW), -1/1 (CW) }
    (A1, A3): 7,
    (A2, A4): 7,
    (A3, A1): 7,
    (A4, A2): 7,

    (A1, OO, A3): 7,
    (A2, OO, A4): 7,
    (A3, OO, A1): 7,
    (A4, OO, A2): 7,

    # Q -- DIA --> Q => { +1/1 (CCW), -1/1 (CW) }
    (Q1, Q3): 8,
    (Q2, Q4): 8,
    (Q3, Q1): 8,
    (Q4, Q2): 8,

    (Q1, OO, Q3): 8,
    (Q2, OO, Q4): 8,
    (Q3, OO, Q1): 8,
    (Q4, OO, Q2): 8,

    # A --- R ---> A => { +1/2 (CCW), -3/2 (CW) }
    (A1, A2): 9,
    (A2, A3): 9,
    (A3, A4): 9,
    (A4, A1): 9,

    (A1, OO, A2): 9,
    (A2, OO, A3): 9,
    (A3, OO, A4): 9,
    (A4, OO, A1): 9,

    # A --- L ---> A => { +3/2 (CCW), -1/2 (CW) }
    (A1, A4): 10,
    (A2, A1): 10,
    (A3, A2): 10,
    (A4, A3): 10,

    (A1, OO, A4): 10,
    (A2, OO, A1): 10,
    (A3, OO, A2): 10,
    (A4, OO, A3): 10,

    # Q --- 1 ---> A => { +3/4 (CCW), -5/4 (CW) }
    (Q1, A3): 11,
    (Q2, A4): 11,
    (Q3, A1): 11,
    (Q4, A2): 11,

    (Q1, OO, A3): 11,
    (Q2, OO, A4): 11,
    (Q3, OO, A1): 11,
    (Q4, OO, A2): 11,

    # Q --- 2 ---> A => { +5/4 (CCW), -3/4 (CW) }
    (Q1, A4): 12,
    (Q2, A1): 12,
    (Q3, A2): 12,
    (Q4, A3): 12,

    (Q1, OO, A4): 12,
    (Q2, OO, A1): 12,
    (Q3, OO, A2): 12,
    (Q4, OO, A3): 12,

    # A --- 1 ---> Q => { +5/4 (CCW), -3/4 (CW) }
    (A1, Q3): 13,
    (A2, Q4): 13,
    (A3, Q1): 13,
    (A4, Q2): 13,

    (A1, OO, Q3): 13,
    (A2, OO, Q4): 13,
    (A3, OO, Q1): 13,
    (A4, OO, Q2): 13,

    # A --- 2 ---> Q => { +3/4 (CCW), -5/4 (CW) }
    (A1, Q2): 14,
    (A2, Q3): 14,
    (A3, Q4): 14,
    (A4, Q1): 14,

    (A1, OO, Q2): 14,
    (A2, OO, Q3): 14,
    (A3, OO, Q4): 14,
    (A4, OO, Q1): 14,

    # Q --> OO --> Q => { +1/2 (CCW), -3/2 (CW) }
    (Q1, OO, Q2): 15,
    (Q2, OO, Q3): 15,
    (Q3, OO, Q4): 15,
    (Q4, OO, Q1): 15,

    # Q --> OO --> Q => { +3/2 (CCW), -1/2 (CW) }
    (Q1, OO, Q4): 16,
    (Q2, OO, Q1): 16,
    (Q3, OO, Q2): 16,
    (Q4, OO, Q3): 16,

    # A --> OO --> A => { +2/1 (CCW), 0 (CW) }
    (A1, OO, A1): 17,
    (A2, OO, A2): 17,
    (A3, OO, A3): 17,
    (A4, OO, A4): 17,

    # Q --> OO --> Q => { +2/1 (CCW), 0 (CW) }
    (Q1, OO, Q1): 18,
    (Q2, OO, Q2): 18,
    (Q3, OO, Q3): 18,
    (Q4, OO, Q4): 18,
}

_values = {
     1: [(+1, 4)],
     2: [(-1, 4)],
     3: [(+1, 4)],
     4: [(-1, 4)],
    -1: [(+9, 4), (+1, 4)],
    -2: [(+7, 4), (-1, 4)],
    -3: [(+9, 4), (+1, 4)],
    -4: [(+7, 4), (-1, 4)],
    +5: [(+1, 2)],
    -5: [(-1, 2)],
     7: [(+1, 1), (-1, 1)],
     8: [(+1, 1), (-1, 1)],
     9: [(+1, 2), (-3, 2)],
    10: [(+3, 2), (-1, 2)],
    11: [(+3, 4), (-5, 4)],
    12: [(+5, 4), (-3, 4)],
    13: [(+5, 4), (-3, 4)],
    14: [(+3, 4), (-5, 4)],
    15: [(+1, 2), (-3, 2)],
    16: [(+3, 2), (-1, 2)],
    17: [(+2, 1), ( 0, 1)],
    18: [(+2, 1), ( 0, 1)],
}

def _classify_point(re, im):
    """Return the half-axis (or origin) on which (re, im) point is located. """
    if not re and not im:
        return OO

    if not re:
        if im > 0:
            return A2
        else:
            return A4
    elif not im:
        if re > 0:
            return A1
        else:
            return A3

def _intervals_to_quadrants(intervals, f1, f2, s, t, F):
    """Generate a sequence of extended quadrants from a list of critical points. """
    if not intervals:
        return []

    Q = []

    if not f1:
        (a, b), _, _ = intervals[0]

        if a == b == s:
            if len(intervals) == 1:
                if dup_eval(f2, t, F) > 0:
                    return [OO, A2]
                else:
                    return [OO, A4]
            else:
                (a, _), _, _ = intervals[1]

                if dup_eval(f2, (s+a)/2, F) > 0:
                    Q.extend([OO, A2])
                    f2_sgn = +1
                else:
                    Q.extend([OO, A4])
                    f2_sgn = -1

                intervals = intervals[1:]
        else:
            if dup_eval(f2, s, F) > 0:
                Q.append(A2)
                f2_sgn = +1
            else:
                Q.append(A4)
                f2_sgn = -1

        for (a, _), indices, _ in intervals:
            Q.append(OO)

            if indices[1] % 2 == 1:
                f2_sgn = -f2_sgn

            if a != t:
                if f2_sgn > 0:
                    Q.append(A2)
                else:
                    Q.append(A4)

        return Q

    if not f2:
        (a, b), _, _ = intervals[0]

        if a == b == s:
            if len(intervals) == 1:
                if dup_eval(f1, t, F) > 0:
                    return [OO, A1]
                else:
                    return [OO, A3]
            else:
                (a, _), _, _ = intervals[1]

                if dup_eval(f1, (s+a)/2, F) > 0:
                    Q.extend([OO, A1])
                    f1_sgn = +1
                else:
                    Q.extend([OO, A3])
                    f1_sgn = -1

                intervals = intervals[1:]
        else:
            if dup_eval(f1, s, F) > 0:
                Q.append(A1)
                f1_sgn = +1
            else:
                Q.append(A3)
                f1_sgn = -1

        for (a, _), indices, _ in intervals:
            Q.append(OO)

            if indices[0] % 2 == 1:
                f1_sgn = -f1_sgn

            if a != t:
                if f1_sgn > 0:
                    Q.append(A1)
                else:
                    Q.append(A3)

        return Q

    re = dup_eval(f1, s, F)
    im = dup_eval(f2, s, F)

    if not re or not im:
        Q.append(_classify_point(re, im))

        if len(intervals) == 1:
            re = dup_eval(f1, t, F)
            im = dup_eval(f2, t, F)
        else:
            (a, _), _, _ = intervals[1]

            re = dup_eval(f1, (s+a)/2, F)
            im = dup_eval(f2, (s+a)/2, F)

        intervals = intervals[1:]

    if re > 0:
        f1_sgn = +1
    else:
        f1_sgn = -1

    if im > 0:
        f2_sgn = +1
    else:
        f2_sgn = -1

    sgn = {
        (+1, +1): Q1,
        (-1, +1): Q2,
        (-1, -1): Q3,
        (+1, -1): Q4,
    }

    Q.append(sgn[(f1_sgn, f2_sgn)])

    for (a, b), indices, _ in intervals:
        if a == b:
            re = dup_eval(f1, a, F)
            im = dup_eval(f2, a, F)

            cls = _classify_point(re, im)

            if cls is not None:
                Q.append(cls)

        if 0 in indices:
            if indices[0] % 2 == 1:
                f1_sgn = -f1_sgn

        if 1 in indices:
            if indices[1] % 2 == 1:
                f2_sgn = -f2_sgn

        if not (a == b and b == t):
            Q.append(sgn[(f1_sgn, f2_sgn)])

    return Q

def _traverse_quadrants(Q_L1, Q_L2, Q_L3, Q_L4, exclude=None):
    """Transform sequences of quadrants to a sequence of rules. """
    if exclude is True:
        edges = [1, 1, 0, 0]

        corners = {
            (0, 1): 1,
            (1, 2): 1,
            (2, 3): 0,
            (3, 0): 1,
        }
    else:
        edges = [0, 0, 0, 0]

        corners = {
            (0, 1): 0,
            (1, 2): 0,
            (2, 3): 0,
            (3, 0): 0,
        }

    if exclude is not None and exclude is not True:
        exclude = set(exclude)

        for i, edge in enumerate(['S', 'E', 'N', 'W']):
            if edge in exclude:
                edges[i] = 1

        for i, corner in enumerate(['SW', 'SE', 'NE', 'NW']):
            if corner in exclude:
                corners[((i - 1) % 4, i)] = 1

    QQ, rules = [Q_L1, Q_L2, Q_L3, Q_L4], []

    for i, Q in enumerate(QQ):
        if not Q:
            continue

        if Q[-1] == OO:
            Q = Q[:-1]

        if Q[0] == OO:
            j, Q = (i - 1) % 4, Q[1:]
            qq = (QQ[j][-2], OO, Q[0])

            if qq in _rules_ambiguous:
                rules.append((_rules_ambiguous[qq], corners[(j, i)]))
            else:
                raise NotImplementedError("3 element rule (corner): " + str(qq))

        q1, k = Q[0], 1

        while k < len(Q):
            q2, k = Q[k], k+1

            if q2 != OO:
                qq = (q1, q2)

                if qq in _rules_simple:
                    rules.append((_rules_simple[qq], 0))
                elif qq in _rules_ambiguous:
                    rules.append((_rules_ambiguous[qq], edges[i]))
                else:
                    raise NotImplementedError("2 element rule (inside): " + str(qq))
            else:
                qq, k = (q1, q2, Q[k]), k+1

                if qq in _rules_ambiguous:
                    rules.append((_rules_ambiguous[qq], edges[i]))
                else:
                    raise NotImplementedError("3 element rule (edge): " + str(qq))

            q1 = qq[-1]

    return rules

def _reverse_intervals(intervals):
    """Reverse intervals for traversal from right to left and from top to bottom. """
    return [ ((b, a), indices, f) for (a, b), indices, f in reversed(intervals) ]

def _winding_number(T, field):
    """Compute the winding number of the input polynomial, i.e. the number of roots. """
    return int(sum([ field(*_values[t][i]) for t, i in T ]) / 2)

def dup_count_complex_roots(f, (u, v), (s, t), K, exclude=None):
    """Count all roots in [u + v*I, s + t*I] rectangle using Collins-Krandick algorithm. """
    if not K.is_ZZ and not K.is_QQ:
        raise DomainError("complex root counting is not supported over %s" % K)

    if K.is_ZZ:
        R, F = K, K.get_field()
    else:
        R, F = K.get_ring(), K

    f1, f2 = dup_real_imag(f, K)

    f1 = dmp_convert(f1, 1, K, F)
    f2 = dmp_convert(f2, 1, K, F)

    f1L1F = dmp_eval_in(f1, v, 1, 1, F)
    f2L1F = dmp_eval_in(f2, v, 1, 1, F)

    _, f1L1R = dup_clear_denoms(f1L1F, F, R, convert=True)
    _, f2L1R = dup_clear_denoms(f2L1F, F, R, convert=True)

    f1L2F = dmp_eval_in(f1, s, 0, 1, F)
    f2L2F = dmp_eval_in(f2, s, 0, 1, F)

    _, f1L2R = dup_clear_denoms(f1L2F, F, R, convert=True)
    _, f2L2R = dup_clear_denoms(f2L2F, F, R, convert=True)

    f1L3F = dmp_eval_in(f1, t, 1, 1, F)
    f2L3F = dmp_eval_in(f2, t, 1, 1, F)

    _, f1L3R = dup_clear_denoms(f1L3F, F, R, convert=True)
    _, f2L3R = dup_clear_denoms(f2L3F, F, R, convert=True)

    f1L4F = dmp_eval_in(f1, u, 0, 1, F)
    f2L4F = dmp_eval_in(f2, u, 0, 1, F)

    _, f1L4R = dup_clear_denoms(f1L4F, F, R, convert=True)
    _, f2L4R = dup_clear_denoms(f2L4F, F, R, convert=True)

    S_L1 = [f1L1R, f2L1R]
    S_L2 = [f1L2R, f2L2R]
    S_L3 = [f1L3R, f2L3R]
    S_L4 = [f1L4R, f2L4R]

    I_L1 = dup_isolate_real_roots_list(S_L1, R, inf=u, sup=s, fast=True, basis=True, strict=True)
    I_L2 = dup_isolate_real_roots_list(S_L2, R, inf=v, sup=t, fast=True, basis=True, strict=True)
    I_L3 = dup_isolate_real_roots_list(S_L3, R, inf=u, sup=s, fast=True, basis=True, strict=True)
    I_L4 = dup_isolate_real_roots_list(S_L4, R, inf=v, sup=t, fast=True, basis=True, strict=True)

    I_L3 = _reverse_intervals(I_L3)
    I_L4 = _reverse_intervals(I_L4)

    Q_L1 = _intervals_to_quadrants(I_L1, f1L1F, f2L1F, u, s, F)
    Q_L2 = _intervals_to_quadrants(I_L2, f1L2F, f2L2F, v, t, F)
    Q_L3 = _intervals_to_quadrants(I_L3, f1L3F, f2L3F, s, u, F)
    Q_L4 = _intervals_to_quadrants(I_L4, f1L4F, f2L4F, t, v, F)

    T = _traverse_quadrants(Q_L1, Q_L2, Q_L3, Q_L4, exclude=exclude)

    return _winding_number(T, F)

def _find_smallest_rectangle(rectangles):
    """Find a rectangle of minimum area for bisection. """
    min_area, j = None, None

    for i, (_, (u, v), (s, t), _, _, _, _) in enumerate(rectangles):
        area = (s - u)*(t - v)

        if min_area is None or area < min_area:
            min_area, j = area, i

    return j

def _rectangle_small_p((u, v), (s, t), eps):
    """Return ``True`` if the given rectangle is small enough. """
    if eps is not None:
        return s - u < eps and t - v < eps
    else:
        return True

def dup_isolate_complex_roots_sqf(f, K, eps=None, inf=None, sup=None):
    """Isolate complex roots of a square-free polynomial using Collins-Krandick algorithm. """
    if not K.is_ZZ and not K.is_QQ:
        raise DomainError("isolation of complex roots is not supported over %s" % K)

    if dup_degree(f) <= 0:
        return []

    if K.is_ZZ:
        R, F = K, K.get_field()
    else:
        R, F = K.get_ring(), K

    f = dup_convert(f, K, F)

    n, lc = dup_degree(f), abs(dup_LC(f, F))
    B = 2*max([ F.quo(abs(c), lc) for c in f ])

    (u, v), (s, t) = (-B, F.zero), (B, B)

    if inf is not None:
        u = inf

    if sup is not None:
        s = sup

    if v < 0 or t <= v or s <= u:
        raise ValueError("not a valid complex isolation rectangle")

    f1, f2 = dup_real_imag(f, F)

    f1L1F = dmp_eval_in(f1, v, 1, 1, F)
    f2L1F = dmp_eval_in(f2, v, 1, 1, F)

    _, f1L1R = dup_clear_denoms(f1L1F, F, R, convert=True)
    _, f2L1R = dup_clear_denoms(f2L1F, F, R, convert=True)

    f1L2F = dmp_eval_in(f1, s, 0, 1, F)
    f2L2F = dmp_eval_in(f2, s, 0, 1, F)

    _, f1L2R = dup_clear_denoms(f1L2F, F, R, convert=True)
    _, f2L2R = dup_clear_denoms(f2L2F, F, R, convert=True)

    f1L3F = dmp_eval_in(f1, t, 1, 1, F)
    f2L3F = dmp_eval_in(f2, t, 1, 1, F)

    _, f1L3R = dup_clear_denoms(f1L3F, F, R, convert=True)
    _, f2L3R = dup_clear_denoms(f2L3F, F, R, convert=True)

    f1L4F = dmp_eval_in(f1, u, 0, 1, F)
    f2L4F = dmp_eval_in(f2, u, 0, 1, F)

    _, f1L4R = dup_clear_denoms(f1L4F, F, R, convert=True)
    _, f2L4R = dup_clear_denoms(f2L4F, F, R, convert=True)

    S_L1 = [f1L1R, f2L1R]
    S_L2 = [f1L2R, f2L2R]
    S_L3 = [f1L3R, f2L3R]
    S_L4 = [f1L4R, f2L4R]

    I_L1 = dup_isolate_real_roots_list(S_L1, R, inf=u, sup=s, fast=True, strict=True, basis=True)
    I_L2 = dup_isolate_real_roots_list(S_L2, R, inf=v, sup=t, fast=True, strict=True, basis=True)
    I_L3 = dup_isolate_real_roots_list(S_L3, R, inf=u, sup=s, fast=True, strict=True, basis=True)
    I_L4 = dup_isolate_real_roots_list(S_L4, R, inf=v, sup=t, fast=True, strict=True, basis=True)

    I_L3 = _reverse_intervals(I_L3)
    I_L4 = _reverse_intervals(I_L4)

    Q_L1 = _intervals_to_quadrants(I_L1, f1L1F, f2L1F, u, s, F)
    Q_L2 = _intervals_to_quadrants(I_L2, f1L2F, f2L2F, v, t, F)
    Q_L3 = _intervals_to_quadrants(I_L3, f1L3F, f2L3F, s, u, F)
    Q_L4 = _intervals_to_quadrants(I_L4, f1L4F, f2L4F, t, v, F)

    T = _traverse_quadrants(Q_L1, Q_L2, Q_L3, Q_L4)
    N = _winding_number(T, F)

    if not N:
        return []

    I = (I_L1, I_L2, I_L3, I_L4)
    Q = (Q_L1, Q_L2, Q_L3, Q_L4)

    F1 = (f1L1F, f1L2F, f1L3F, f1L4F)
    F2 = (f2L1F, f2L2F, f2L3F, f2L4F)

    rectangles, roots = [(N, (u, v), (s, t), I, Q, F1, F2)], []

    while rectangles:
        i = _find_smallest_rectangle(rectangles)
        N, (u, v), (s, t), II, QQ, F1, F2 = rectangles.pop(i)

        I_L1, I_L2, I_L3, I_L4 = II
        Q_L1, Q_L2, Q_L3, Q_L4 = QQ

        f1L1F, f1L2F, f1L3F, f1L4F = F1
        f2L1F, f2L2F, f2L3F, f2L4F = F2

        if s - u > t - v:
            x = (u + s) / 2

            f1VF = dmp_eval_in(f1, x, 0, 1, F)
            f2VF = dmp_eval_in(f2, x, 0, 1, F)

            _, f1VR = dup_clear_denoms(f1VF, F, R, convert=True)
            _, f2VR = dup_clear_denoms(f2VF, F, R, convert=True)

            I_V = dup_isolate_real_roots_list([f1VR, f2VR], R, inf=v, sup=t, fast=True, strict=True, basis=True)

            I_L1_L, I_L1_R = [], []
            I_L2_L, I_L2_R = I_V, I_L2
            I_L3_L, I_L3_R = [], []
            I_L4_L, I_L4_R = I_L4, _reverse_intervals(I_V)

            for I in I_L1:
                (a, b), indices, h = I

                if a == b:
                    if a == x:
                        I_L1_L.append(I)
                        I_L1_R.append(I)
                    elif a < x:
                        I_L1_L.append(I)
                    else:
                        I_L1_R.append(I)
                else:
                    if b <= x:
                        I_L1_L.append(I)
                    elif a >= x:
                        I_L1_R.append(I)
                    else:
                        a, b = dup_refine_real_root(h, a, b, R, disjoint=x, fast=True)

                        if b <= x:
                            I_L1_L.append(((a, b), indices, h))
                        else:
                            I_L1_R.append(((a, b), indices, h))

            for I in I_L3:
                (b, a), indices, h = I

                if a == b:
                    if a == x:
                        I_L3_L.append(I)
                        I_L3_R.append(I)
                    elif a < x:
                        I_L3_L.append(I)
                    else:
                        I_L3_R.append(I)
                else:
                    if b <= x:
                        I_L3_L.append(I)
                    elif a >= x:
                        I_L3_R.append(I)
                    else:
                        a, b = dup_refine_real_root(h, a, b, R, disjoint=x, fast=True)

                        if b <= x:
                            I_L3_L.append(((b, a), indices, h))
                        else:
                            I_L3_R.append(((b, a), indices, h))

            Q_L1_L = _intervals_to_quadrants(I_L1_L, f1L1F, f2L1F, u, x, F)
            Q_L2_L = _intervals_to_quadrants(I_L2_L, f1VF,  f2VF,  v, t, F)
            Q_L3_L = _intervals_to_quadrants(I_L3_L, f1L3F, f2L3F, x, u, F)
            Q_L4_L = Q_L4

            Q_L1_R = _intervals_to_quadrants(I_L1_R, f1L1F, f2L1F, x, s, F)
            Q_L2_R = Q_L2
            Q_L3_R = _intervals_to_quadrants(I_L3_R, f1L3F, f2L3F, s, x, F)
            Q_L4_R = _intervals_to_quadrants(I_L4_R, f1VF,  f2VF,  t, v, F)

            T_L = _traverse_quadrants(Q_L1_L, Q_L2_L, Q_L3_L, Q_L4_L, exclude=True)
            T_R = _traverse_quadrants(Q_L1_R, Q_L2_R, Q_L3_R, Q_L4_R, exclude=True)

            N_L = _winding_number(T_L, F)
            N_R = _winding_number(T_R, F)

            I_L = (I_L1_L, I_L2_L, I_L3_L, I_L4_L)
            Q_L = (Q_L1_L, Q_L2_L, Q_L3_L, Q_L4_L)

            I_R = (I_L1_R, I_L2_R, I_L3_R, I_L4_R)
            Q_R = (Q_L1_R, Q_L2_R, Q_L3_R, Q_L4_R)

            F1_L = (f1L1F, f1VF, f1L3F, f1L4F)
            F2_L = (f2L1F, f2VF, f2L3F, f2L4F)

            F1_R = (f1L1F, f1L2F, f1L3F, f1VF)
            F2_R = (f2L1F, f2L2F, f2L3F, f2VF)

            if N_L >= 1:
                a, b = (u, v), (x, t)

                if N_L == 1 and _rectangle_small_p(a, b, eps):
                    roots.append((a, b))
                else:
                    rectangles.append((N_L, a, b, I_L, Q_L, F1_L, F2_L))

            if N_R >= 1:
                a, b = (x, v), (s, t)

                if N_R == 1 and _rectangle_small_p(a, b, eps):
                    roots.append((a, b))
                else:
                    rectangles.append((N_R, a, b, I_R, Q_R, F1_R, F2_R))
        else:
            y = (v + t) / 2

            f1HF = dmp_eval_in(f1, y, 1, 1, F)
            f2HF = dmp_eval_in(f2, y, 1, 1, F)

            _, f1HR = dup_clear_denoms(f1HF, F, R, convert=True)
            _, f2HR = dup_clear_denoms(f2HF, F, R, convert=True)

            I_H = dup_isolate_real_roots_list([f1HR, f2HR], R, inf=u, sup=s, fast=True, strict=True, basis=True)

            I_L1_B, I_L1_U = I_L1, I_H
            I_L2_B, I_L2_U = [], []
            I_L3_B, I_L3_U = _reverse_intervals(I_H), I_L3
            I_L4_B, I_L4_U = [], []

            for I in I_L2:
                (a, b), indices, h = I

                if a == b:
                    if a == y:
                        I_L2_B.append(I)
                        I_L2_U.append(I)
                    elif a < y:
                        I_L2_B.append(I)
                    else:
                        I_L2_U.append(I)
                else:
                    if b <= y:
                        I_L2_B.append(I)
                    elif a >= y:
                        I_L2_U.append(I)
                    else:
                        a, b = dup_refine_real_root(h, a, b, R, disjoint=y, fast=True)

                        if b <= y:
                            I_L2_B.append(((a, b), indices, h))
                        else:
                            I_L2_U.append(((a, b), indices, h))

            for I in I_L4:
                (b, a), indices, h = I

                if a == b:
                    if a == y:
                        I_L4_B.append(I)
                        I_L4_U.append(I)
                    elif a < y:
                        I_L4_B.append(I)
                    else:
                        I_L4_U.append(I)
                else:
                    if b <= y:
                        I_L4_B.append(I)
                    elif a >= y:
                        I_L4_U.append(I)
                    else:
                        a, b = dup_refine_real_root(h, a, b, R, disjoint=y, fast=True)

                        if b <= y:
                            I_L4_B.append(((b, a), indices, h))
                        else:
                            I_L4_U.append(((b, a), indices, h))

            Q_L1_B = Q_L1
            Q_L2_B = _intervals_to_quadrants(I_L2_B, f1L2F, f2L2F, v, y, F)
            Q_L3_B = _intervals_to_quadrants(I_L3_B, f1HF,  f2HF,  s, u, F)
            Q_L4_B = _intervals_to_quadrants(I_L4_B, f1L4F, f2L4F, y, v, F)

            Q_L1_U = _intervals_to_quadrants(I_L1_U, f1HF,  f2HF,  u, s, F)
            Q_L2_U = _intervals_to_quadrants(I_L2_U, f1L2F, f2L2F, y, t, F)
            Q_L3_U = Q_L3
            Q_L4_U = _intervals_to_quadrants(I_L4_U, f1L4F, f2L4F, t, y, F)

            T_B = _traverse_quadrants(Q_L1_B, Q_L2_B, Q_L3_B, Q_L4_B, exclude=True)
            T_U = _traverse_quadrants(Q_L1_U, Q_L2_U, Q_L3_U, Q_L4_U, exclude=True)

            N_B = _winding_number(T_B, F)
            N_U = _winding_number(T_U, F)

            I_B = (I_L1_B, I_L2_B, I_L3_B, I_L4_B)
            Q_B = (Q_L1_B, Q_L2_B, Q_L3_B, Q_L4_B)

            I_U = (I_L1_U, I_L2_U, I_L3_U, I_L4_U)
            Q_U = (Q_L1_U, Q_L2_U, Q_L3_U, Q_L4_U)

            F1_B = (f1L1F, f1L2F, f1HF, f1L4F)
            F2_B = (f2L1F, f2L2F, f2HF, f2L4F)

            F1_U = (f1HF, f1L2F, f1L3F, f1L4F)
            F2_U = (f2HF, f2L2F, f2L3F, f2L4F)

            if N_B >= 1:
                a, b = (u, v), (s, y)

                if N_B == 1 and _rectangle_small_p(a, b, eps):
                    roots.append((a, b))
                else:
                    rectangles.append((N_B, a, b, I_B, Q_B, F1_B, F2_B))

            if N_U >= 1:
                a, b = (u, y), (s, t)

                if N_U == 1 and _rectangle_small_p(a, b, eps):
                    roots.append((a, b))
                else:
                    rectangles.append((N_U, a, b, I_U, Q_U, F1_U, F2_U))

    _roots, roots = sorted(roots, key=operator.itemgetter(0)), []

    for (u, v), (s, t) in _roots:
        roots.extend([((u, v), (s, t)), ((u, -t), (s, -v))])

    return roots

def dup_isolate_all_roots_sqf(f, K, eps=None, inf=None, sup=None, fast=False):
    """Isolate real and complex roots of a square-free polynomial ``f``. """
    return (dup_isolate_real_roots_sqf(f, K, eps=eps, inf=inf, sup=sup, fast=fast),
            dup_isolate_complex_roots_sqf(f, K, eps=eps, inf=inf, sup=sup))

def dup_isolate_all_roots(f, K, eps=None, inf=None, sup=None, fast=False):
    """Isolate real and complex roots of a non-square-free polynomial ``f``. """
    if not K.is_ZZ and not K.is_QQ:
        raise DomainError("isolation of real and complex roots is not supported over %s" % K)

    _, factors = dup_sqf_list(f, K)

    if len(factors) == 1:
        ((f, k),) = factors

        real_part, complex_part = dup_isolate_all_roots_sqf(f, K, eps=eps, inf=inf, sup=sup, fast=fast)

        real_part = [ ((a, b), k) for (a, b) in real_part ]
        complex_part = [ ((a, b), k) for (a, b) in complex_part ]

        return real_part, complex_part
    else:
        raise NotImplementedError("only trivial square-free polynomials are supported")

