
from sympy.core.numbers import igcdex

def crt(m, v, symmetric=False):
    """Chinese Remainder Theorem.

       The integers in m are assumed to be pairwise coprime.  The output
       is then an integer f, such that f = v_i mod m_i for each pair out
       of v and m.

    """
    mm = 1

    for m_i in m:
        mm *= m_i

    result = 0

    for m_i, v_i in zip(m, v):
        e = mm // m_i
        s, t, g = igcdex(e, m_i)
        c = v_i*s % m_i
        result += c*e

    result %= mm

    if symmetric:
        if result <= mm//2:
            return result
        else:
            return result - mm
    else:
        return result

def crt1(m):
    """First part of chines remainder theorem, for multiple application. """
    mm, e, s = 1, [], []

    for m_i in m:
        mm *= m_i

    for m_i in m:
        e.append(mm // m_i)
        s.append(igcdex(e[-1], m_i)[0])

    return mm, e, s

def crt2(m, v, mm, e, s, symmetric=False):
    """Second part of Chinese Remainder Theorem, for multiple application. """
    result = 0

    for m_i, v_i, e_i, s_i in zip(m, v, e, s):
        c = v_i*s_i % m_i
        result += c*e_i

    result %= mm

    if symmetric:
        if result <= mm // 2:
            return result
        else:
            return result - mm
    else:
        return result
