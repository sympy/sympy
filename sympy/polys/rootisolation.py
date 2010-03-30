"""Tools for real and complex root isolation and refinement. """

from sympy.polys.densebasic import (
    dup_degree, dmp_convert, dup_terms_gcd,
)

from sympy.polys.densetools import (
    dup_isolate_real_roots_list,
    dup_eval, dmp_eval_in,
    dup_clear_denoms,
    dup_real_imag,
    dup_inner_gcd,
)

from sympy.polys.polyerrors import (
    DomainError,
)

OO = 'OO' # Origin of (re, im) coordinate system

Q1 = 'Q1' # Quadrant #1 (++): re > 0 and im > 0
Q2 = 'Q2' # Quadrant #2 (-+): re < 0 and im > 0
Q3 = 'Q3' # Quadrant #3 (--): re < 0 and im < 0
Q4 = 'Q4' # Quadrant #4 (+-): re > 0 and im < 0

A1 = 'A1' # Axis #1 (+0): re > 0 and im = 0
A2 = 'A2' # Axis #2 (0+): re = 0 and im > 0
A3 = 'A3' # Axis #3 (-0): re < 0 and im = 0
A4 = 'A4' # Axis #4 (0-): re = 0 and im < 0

_rule_values = {
     1: (+1, 4),
     2: (-1, 4),
   200: (+7, 4),
     3: (+1, 4),
     4: (-1, 4),
   400: (+7, 4),
     5: (+1, 2),
     6: (-1, 2),
     7: (+1, 1),
     8: (+1, 1),
     9: (+1, 2),
    10: (+3, 2),
    11: (+3, 4),
    12: (+5, 4),
    13: (+5, 4),
    14: (+3, 4),
    15: (+1, 2),
    16: (+3, 2),
    17: (+2, 1),
    18: (+2, 1),
}

_quadrant_rules = {
    # A -- CCW --> Q => +1/4, +9/4
    (A1, Q1): 1,
    (A2, Q2): 1,
    (A3, Q3): 1,
    (A4, Q4): 1,

    # A --  CW --> Q => +7/4, -1,4
    (A1, Q4): 2,
    (A2, Q1): 2,
    (A3, Q2): 2,
    (A4, Q3): 2,

    (A1, OO, Q4): 200,
    (A2, OO, Q1): 200,
    (A3, OO, Q2): 200,
    (A4, OO, Q3): 200,

    # Q -- CCW --> A => +9/4, +1/4
    (Q1, A2): 3,
    (Q2, A3): 3,
    (Q3, A4): 3,
    (Q4, A1): 3,

    # Q --  CW --> A => +7/4, -1/4
    (Q1, A1): 4,
    (Q2, A2): 4,
    (Q3, A3): 4,
    (Q4, A4): 4,

    (Q1, OO, A1): 400,
    (Q2, OO, A2): 400,
    (Q3, OO, A3): 400,
    (Q4, OO, A4): 400,

    # Q -- CCW --> Q => +1/2
    (Q1, Q2): 5,
    (Q2, Q3): 5,
    (Q3, Q4): 5,
    (Q4, Q1): 5,

    # Q --  CW --> Q => -1/2
    (Q1, Q4): 6,
    (Q2, Q1): 6,
    (Q3, Q2): 6,
    (Q4, Q3): 6,

    # A ---------> A => +1, -1
    (A1, A3): 7,
    (A2, A4): 7,
    (A3, A1): 7,
    (A4, A2): 7,

    (A1, OO, A3): 7,
    (A2, OO, A4): 7,
    (A3, OO, A1): 7,
    (A4, OO, A2): 7,

    # Q -- DIA --> Q => +1, -1
    (Q1, Q3): 8,
    (Q2, Q4): 8,
    (Q3, Q1): 8,
    (Q4, Q2): 8,

    (Q1, OO, Q3): 8,
    (Q2, OO, Q4): 8,
    (Q3, OO, Q1): 8,
    (Q4, OO, Q2): 8,

    # A --- R ---> A => +1/2, -3/2
    (A1, A2): 9,
    (A2, A3): 9,
    (A3, A4): 9,
    (A4, A1): 9,

    (A1, OO, A2): 9,
    (A2, OO, A3): 9,
    (A3, OO, A4): 9,
    (A4, OO, A1): 9,

    # A --- L ---> A => -1/2, +3/2
    (A1, A4): 10,
    (A2, A1): 10,
    (A3, A2): 10,
    (A4, A3): 10,

    (A1, OO, A4): 10,
    (A2, OO, A1): 10,
    (A3, OO, A2): 10,
    (A4, OO, A3): 10,

    # Q --- 1 ---> A => +3/4, -5/4
    (Q1, A3): 11,
    (Q2, A4): 11,
    (Q3, A1): 11,
    (Q4, A2): 11,

    (Q1, OO, A3): 11,
    (Q2, OO, A4): 11,
    (Q3, OO, A1): 11,
    (Q4, OO, A2): 11,

    # Q --- 2 ---> A => +5/4, -3/4
    (Q1, A4): 12,
    (Q2, A1): 12,
    (Q3, A2): 12,
    (Q4, A3): 12,

    # A --- 1 ---> Q => +5/4, -3/4
    (A1, Q3): 13,
    (A2, Q4): 13,
    (A3, Q1): 13,
    (A4, Q2): 13,

    # A --- 2 ---> Q => +3/4, -5/4
    (A1, Q2): 14,
    (A2, Q3): 14,
    (A3, Q4): 14,
    (A4, Q1): 14,

    # Q --> OO --> Q => +1/2 (CCW)
    (Q1, OO, Q2): 15,
    (Q2, OO, Q3): 15,
    (Q3, OO, Q4): 15,
    (Q4, OO, Q1): 15,

    # Q --> OO --> Q => -3/2 (CW?)
    (Q1, OO, Q4): 16,
    (Q2, OO, Q1): 16,
    (Q3, OO, Q2): 16,
    (Q4, OO, Q3): 16,

    # A --> OO --> A => +2 (SELF)
    (A1, OO, A1): 17,
    (A2, OO, A2): 17,
    (A3, OO, A3): 17,
    (A4, OO, A4): 17,

    # Q --> OO --> Q => +2 (SELF)
    (Q1, OO, Q1): 18,
    (Q2, OO, Q2): 18,
    (Q3, OO, Q3): 18,
    (Q4, OO, Q4): 18,
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

def _intervals_to_quadrands(intervals, f1, f2, s, t, F):
    """Generate a sequence of extended quadrants from a list of critical points. """
    if not intervals:
        return []

    Q = []

    if not f1:
        (a, _), _ = intervals[0]

        if a == s:
            if len(intervals) == 1:
                if dup_eval(f2, t, F) > 0:
                    return [OO, A2]
                else:
                    return [OO, A4]
            else:
                (a, _), _ = intervals[1]

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

        for (a, _), indices in intervals:
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
        (a, _), _ = intervals[0]

        if a == s:
            if len(intervals) == 1:
                if dup_eval(f1, t, F) > 0:
                    return [OO, A1]
                else:
                    return [OO, A3]
            else:
                (a, _), _ = intervals[1]

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

        for (a, _), indices in intervals:
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
            (a, _), _ = intervals[1]

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

    for (a, b), indices in intervals:
        if a == b:
            re = dup_eval(f1, a, F)
            im = dup_eval(f2, a, F)

            Q.append(_classify_point(re, im))

        if 0 in indices:
            if indices[0] % 2 == 1:
                f1_sgn = -f1_sgn

        if 1 in indices:
            if indices[1] % 2 == 1:
                f2_sgn = -f2_sgn

        if not (a == b and b == t):
            Q.append(sgn[(f1_sgn, f2_sgn)])

    return Q

def _traverse_quadrants(Q):
    """Transform a sequence of quadrants to a sequence of rules. """
    if not Q:
        return []

    while Q[0] == OO:
        Q = Q[1:] + [Q[0]]

    Q_new = [Q[0]]

    for q2 in Q[1:]:
        q1 = Q_new[-1]

        if q1 != q2:
            Q_new.append(q2)

    Q = Q_new

    if Q[-1] == OO:
        Q.append(Q[0])

    q1, i, rules = Q[0], 1, []

    while i < len(Q):
        q2 = Q[i]

        if i+1 < len(Q):
            q3 = Q[i+1]
            qq = q1, q2, q3

            if qq in _quadrant_rules:
                rules.append(_quadrant_rules[qq])
                q1 = q3
                i += 2
                continue

        qq = q1, q2

        if qq in _quadrant_rules:
            rules.append(_quadrant_rules[qq])
            q1 = q2
            i += 1
            continue

        raise RuntimeError("no such quadrant rule: %s" % (qq,))

    return rules

def _winding_number(Q, field):
    """Compute the number of times ``f`` winds around the origin traversing the boundary. """
    return sum([ field(*_rule_values[t]) for t in _traverse_quadrants(Q) ])/2

def dup_count_complex_roots(f, (u, v), (s, t), K):
    """Count all roots in [u + v*I, s + t*I] rectangle using Collins-Krandick algorithm. """
    if not K.is_ZZ and not K.is_QQ:
        raise DomainError("isolation of complex roots not supported over %s" % K)

    if K.is_ZZ:
        R, F = K, K.get_field()
    else:
        R, F = K.get_ring(), K

    k_zero, f = dup_terms_gcd(f, K)

    if u > 0 or s < 0 or v > 0 or t < 0:
        k_zero = 0

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

    I_L1 = dup_isolate_real_roots_list(S_L1, R, inf=u, sup=s, fast=True, strict=True)
    I_L2 = dup_isolate_real_roots_list(S_L2, R, inf=v, sup=t, fast=True, strict=True)
    I_L3 = dup_isolate_real_roots_list(S_L3, R, inf=u, sup=s, fast=True, strict=True)
    I_L4 = dup_isolate_real_roots_list(S_L4, R, inf=v, sup=t, fast=True, strict=True)

    I_L3 = list(reversed(I_L3))
    I_L4 = list(reversed(I_L4))

    Q_L1 = _intervals_to_quadrands(I_L1, f1L1F, f2L1F, u, s, F)
    Q_L2 = _intervals_to_quadrands(I_L2, f1L2F, f2L2F, v, t, F)
    Q_L3 = _intervals_to_quadrands(I_L3, f1L3F, f2L3F, s, u, F)
    Q_L4 = _intervals_to_quadrands(I_L4, f1L4F, f2L4F, t, v, F)

    N = _winding_number(Q_L1 + Q_L2 + Q_L3 + Q_L4, F)

    return k_zero + N
