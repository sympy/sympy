# -*- coding: utf-8 -*-
r"""
Wigner, Clebsch-Gordan, Racah, and Gaunt coefficients

Collection of functions for calculating Wigner 3j, 6j, 9j,
Clebsch-Gordan, Racah as well as Gaunt coefficients exactly, all
evaluating to a rational number times the square root of a rational
number [Rasch03]_.

Please see the description of the individual functions for further
details and examples.

References
==========

.. [Regge58] 'Symmetry Properties of Clebsch-Gordan Coefficients',
  T. Regge, Nuovo Cimento, Volume 10, pp. 544 (1958)
.. [Regge59] 'Symmetry Properties of Racah Coefficients',
  T. Regge, Nuovo Cimento, Volume 11, pp. 116 (1959)
.. [Edmonds74] A. R. Edmonds. Angular momentum in quantum mechanics.
  Investigations in physics, 4.; Investigations in physics, no. 4.
  Princeton, N.J., Princeton University Press, 1957.
.. [Rasch03] J. Rasch and A. C. H. Yu, 'Efficient Storage Scheme for
  Pre-calculated Wigner 3j, 6j and Gaunt Coefficients', SIAM
  J. Sci. Comput. Volume 25, Issue 4, pp. 1416-1428 (2003)
.. [Liberatodebrito82] 'FORTRAN program for the integral of three
  spherical harmonics', A. Liberato de Brito,
  Comput. Phys. Commun., Volume 25, pp. 81-85 (1982)
.. [Homeier96] 'Some Properties of the Coupling Coefficients of Real
  Spherical Harmonics and Their Relation to Gaunt Coefficients',
  H. H. H. Homeier and E. O. Steinborn J. Mol. Struct., Volume 368,
  pp. 31-37 (1996)
.. [Wei99] L. Wei, 'Unified approach for exact calculation of
  angular momentum coupling and recoupling coefficients',
  Computer Physics Communications 120, 222 (1999).
.. [Varshalovich88] D. A. Varshalovich, A. N. Moskalev and
  V. K. Khersonskii, 'Quantum Theory of Angular Momentum',
  (World Scientific, 1988).

Credits and Copyright
=====================

The origin code was taken from Sage with the permission of all authors:

https://groups.google.com/forum/#!topic/sage-devel/M4NZdu-7O38

Now it has been rewritten using new formulations base on [Wei99]_ and
[Varshalovich88]_.

Authors
=======

- Jens Rasch (2009-03-24): initial version for Sage

- Jens Rasch (2009-05-31): updated to sage-4.0

- Oscar Gerardo Lazo Arjona (2017-06-18): added Wigner D matrices

- Phil Adam LeMaitre (2022-09-19): added real Gaunt coefficient

- Shaoliang Jin (2025-04-26): rewrote Wigner 3nj and Gaunt coefficient

Copyright (C) 2008 Jens Rasch <jyr2000@gmail.com>

"""
from sympy.concrete.summations import Sum
from sympy.core.add import Add
from sympy.core.function import Function
from sympy.core.numbers import (I, Integer, pi)
from sympy.core.singleton import S
from sympy.core.symbol import Dummy
from sympy.core.sympify import sympify
from sympy.functions.combinatorial.factorials import binomial
from sympy.functions.elementary.complexes import re
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import (cos, sin)
from sympy.functions.special.spherical_harmonics import Ynm
from sympy.matrices.dense import zeros
from sympy.matrices.immutable import ImmutableMatrix
from sympy.utilities.misc import as_int
from math import comb


def _as_int(value):
    """use non-strict conversion to int to handle floats"""
    return as_int(value, strict=False)

def _doubled_int(value):
    """return double of value so that it is always an Python int"""
    try:
        return _as_int(value * 2)
    except ValueError:
        raise ValueError("expecting integer or half-integer, got %s" % value) from None


def _check_dj_couple(dj1, dj2, dj3):
    """check if three angular momenta couple"""
    if (dj1 + dj2 + dj3) % 2 != 0:
        return False
    if dj1 > dj2 + dj3:
        return False
    if dj2 > dj1 + dj3:
        return False
    if dj3 > dj1 + dj2:
        return False
    return True


def wigner_3j(j_1, j_2, j_3, m_1, m_2, m_3, prec=None):
    r"""
    Calculate the Wigner 3j symbol `\operatorname{Wigner3j}(j_1,j_2,j_3,m_1,m_2,m_3)`.

    Parameters
    ==========

    j_1, j_2, j_3, m_1, m_2, m_3 :
        Integer or half integer.
    prec :
        Precision, default: ``None``.

    Returns
    =======

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    Examples
    ========

    >>> from sympy.physics.wigner import wigner_3j
    >>> wigner_3j(2, 6, 4, 0, 0, 0)
    sqrt(715)/143
    >>> wigner_3j(2, 6, 4, 0, 0, 1)
    0

    It is an error to have arguments that are not integer or half
    integer values::

        >>> wigner_3j(2.1, 6, 4, 0, 0, 0)
        Traceback (most recent call last):
        ...
        ValueError: expecting integer or half-integer, got 2.1
        >>> wigner_3j(2, 6, 4, 1, 0, -1.1)
        Traceback (most recent call last):
        ...
        ValueError: expecting integer or half-integer, got -1.1

    Notes
    =====

    The Wigner 3j symbol obeys the following symmetry rules:

    - invariant under any permutation of the columns (with the
      exception of a sign change where `J:=j_1+j_2+j_3`):

      .. math::

         \begin{aligned}
         \operatorname{Wigner3j}(j_1,j_2,j_3,m_1,m_2,m_3)
          &=\operatorname{Wigner3j}(j_3,j_1,j_2,m_3,m_1,m_2) \\
          &=\operatorname{Wigner3j}(j_2,j_3,j_1,m_2,m_3,m_1) \\
          &=(-1)^J \operatorname{Wigner3j}(j_3,j_2,j_1,m_3,m_2,m_1) \\
          &=(-1)^J \operatorname{Wigner3j}(j_1,j_3,j_2,m_1,m_3,m_2) \\
          &=(-1)^J \operatorname{Wigner3j}(j_2,j_1,j_3,m_2,m_1,m_3)
         \end{aligned}

    - invariant under space inflection, i.e.

      .. math::

         \operatorname{Wigner3j}(j_1,j_2,j_3,m_1,m_2,m_3)
         =(-1)^J \operatorname{Wigner3j}(j_1,j_2,j_3,-m_1,-m_2,-m_3)

    - symmetric with respect to the 72 additional symmetries based on
      the work by [Regge58]_

    - zero for `j_1`, `j_2`, `j_3` not fulfilling triangle relation

    - zero for `m_1 + m_2 + m_3 \neq 0`

    - zero for violating any one of the conditions

      .. math::

         m_1  \in \{-|j_1|, \ldots, |j_1|\},
         m_2  \in \{-|j_2|, \ldots, |j_2|\},
         m_3  \in \{-|j_3|, \ldots, |j_3|\}

    Algorithm
    =========

    This function uses the algorithm of [Wei99]_ to calculate the
    value of the 3j symbol exactly. Note that the formula contains
    alternating sums over large binomials and is therefore unsuitable
    for finite precision arithmetic and only useful for a computer
    algebra system [Rasch03]_.

    Authors
    =======

    - Jens Rasch (2009-03-24): initial version
    """

    dj1, dj2, dj3, dm1, dm2, dm3 = \
        map(_doubled_int, [j_1, j_2, j_3, m_1, m_2, m_3])


    if dm1 + dm2 + dm3 != 0:
        return S.Zero
    if not _check_dj_couple(dj1, dj2, dj3):
        return S.Zero
    if (abs(dm1) > dj1) or (abs(dm2) > dj2) or (abs(dm3) > dj3):
        return S.Zero
    if not ((dj1 - dm1) % 2 == 0 and \
            (dj2 - dm2) % 2 == 0 and \
            (dj3 - dm3) % 2 == 0):
        return S.Zero

    sumj = (dj1 + dj2 + dj3) // 2
    jm1 = sumj - dj1
    jm2 = sumj - dj2
    jm3 = sumj - dj3
    j1mm1 = (dj1 - dm1) // 2
    j2mm2 = (dj2 - dm2) // 2
    j3mm3 = (dj3 - dm3) // 2
    j1pm1 = (dj1 + dm1) // 2

    imin = max(0, j1pm1 - jm2, j2mm2 - jm1)
    imax = min(jm3, j1pm1, j2mm2)
    sumres = 0
    for ii in range(imin, imax + 1):
        ti = comb(jm3, ii) * comb(jm2, j1pm1 - ii) * \
            comb(jm1, j2mm2 - ii)
        sumres = ti - sumres
    if sumres == 0:
        return S.Zero

    res_num = comb(dj1, jm2) * comb(dj2, jm1)
    res_den = comb(sumj, jm3) * comb(dj1, j1mm1) * \
        comb(dj2, j2mm2) * comb(dj3, j3mm3) * (sumj + 1)
    ressqrt = sqrt(Integer(res_num) / res_den)

    phase = (-1) ** (dj1 + (dj3 + dm3) // 2 + imax)
    res = phase * ressqrt * sumres
    if prec:
        res = res.evalf(prec)
    return res


def clebsch_gordan(j_1, j_2, j_3, m_1, m_2, m_3, prec=None):
    r"""
    Calculates the Clebsch-Gordan coefficient.
    `\left\langle j_1 m_1 \; j_2 m_2 | j_3 m_3 \right\rangle`.

    The reference for this function is [Edmonds74]_.

    Parameters
    ==========

    j_1, j_2, j_3, m_1, m_2, m_3 :
        Integer or half integer.
    prec :
        Precision, default: ``None``.

    Returns
    =======

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    Examples
    ========

    >>> from sympy import S
    >>> from sympy.physics.wigner import clebsch_gordan
    >>> clebsch_gordan(S(3)/2, S(1)/2, 2, S(3)/2, S(1)/2, 2)
    1
    >>> clebsch_gordan(S(3)/2, S(1)/2, 1, S(3)/2, -S(1)/2, 1)
    sqrt(3)/2
    >>> clebsch_gordan(S(3)/2, S(1)/2, 1, -S(1)/2, S(1)/2, 0)
    -sqrt(2)/2

    Notes
    =====

    The Clebsch-Gordan coefficient will be evaluated via its relation
    to Wigner 3j symbols:

    .. math::

        \left\langle j_1 m_1 \; j_2 m_2 | j_3 m_3 \right\rangle
        =(-1)^{j_1-j_2+m_3} \sqrt{2j_3+1}
        \operatorname{Wigner3j}(j_1,j_2,j_3,m_1,m_2,-m_3)

    See also the documentation on Wigner 3j symbols which exhibit much
    higher symmetry relations than the Clebsch-Gordan coefficient.

    Authors
    =======

    - Jens Rasch (2009-03-24): initial version
    """
    w = wigner_3j(j_1, j_2, j_3, m_1, m_2, -m_3, prec=prec)
    return (-1) ** (j_1 - j_2 + m_3) * sqrt(Integer(2 * j_3 + 1)) * w


def wigner_6j(j_1, j_2, j_3, j_4, j_5, j_6, prec=None):
    r"""
    Calculate the Wigner 6j symbol `\operatorname{Wigner6j}(j_1,j_2,j_3,j_4,j_5,j_6)`.

    Parameters
    ==========

    j_1, ..., j_6 :
        Integer or half integer.
    prec :
        Precision, default: ``None``.

    Returns
    =======

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    Examples
    ========

    >>> from sympy.physics.wigner import wigner_6j
    >>> wigner_6j(3,3,3,3,3,3)
    -1/14
    >>> wigner_6j(5,5,5,5,5,5)
    1/52

    It is an error to have arguments that are not integer or half
    integer values::

        >>> wigner_6j(0.5,0.5,1.1,0.5,0.5,1.1)
        Traceback (most recent call last):
        ...
        ValueError: expecting integer or half-integer, got 1.1

    Notes
    =====

    The Wigner 6j symbol is related to the Racah symbol but exhibits
    more symmetries as detailed below.

    .. math::

       \operatorname{Wigner6j}(j_1,j_2,j_3,j_4,j_5,j_6)
        =(-1)^{j_1+j_2+j_4+j_5} W(j_1,j_2,j_5,j_4,j_3,j_6)

    The Wigner 6j symbol obeys the following symmetry rules:

    - Wigner 6j symbols are left invariant under any permutation of
      the columns:

      .. math::

         \begin{aligned}
         \operatorname{Wigner6j}(j_1,j_2,j_3,j_4,j_5,j_6)
          &=\operatorname{Wigner6j}(j_3,j_1,j_2,j_6,j_4,j_5) \\
          &=\operatorname{Wigner6j}(j_2,j_3,j_1,j_5,j_6,j_4) \\
          &=\operatorname{Wigner6j}(j_3,j_2,j_1,j_6,j_5,j_4) \\
          &=\operatorname{Wigner6j}(j_1,j_3,j_2,j_4,j_6,j_5) \\
          &=\operatorname{Wigner6j}(j_2,j_1,j_3,j_5,j_4,j_6)
         \end{aligned}

    - They are invariant under the exchange of the upper and lower
      arguments in each of any two columns, i.e.

      .. math::

         \begin{aligned}
         \operatorname{Wigner6j}(j_1,j_2,j_3,j_4,j_5,j_6)
          &=\operatorname{Wigner6j}(j_1,j_5,j_6,j_4,j_2,j_3)\\
          &=\operatorname{Wigner6j}(j_4,j_2,j_6,j_1,j_5,j_3)\\
          &=\operatorname{Wigner6j}(j_4,j_5,j_3,j_1,j_2,j_6)
         \end{aligned}

    - additional 6 symmetries [Regge59]_ giving rise to 144 symmetries
      in total

    - zero if any of the following triples of `j`'s do not fulfill
      the triangle relation

      .. math::

         \{j_1, j_2, j_3\}, \{j_1, j_5, j_6\},
         \{j_4, j_2, j_6\}, \{j_4, j_5, j_3\}

    Algorithm
    =========

    This function uses the algorithm of [Wei99]_ to calculate the
    value of the 6j symbol exactly. Note that the formula contains
    alternating sums over large factorials and is therefore unsuitable
    for finite precision arithmetic and only useful for a computer
    algebra system [Rasch03]_.
    """
    dj1, dj2, dj3, dj4, dj5, dj6 = \
          map(_doubled_int, [j_1, j_2, j_3, j_4, j_5, j_6])
    if not _check_dj_couple(dj1, dj2, dj3):
        return S.Zero
    if not _check_dj_couple(dj1, dj5, dj6):
        return S.Zero
    if not _check_dj_couple(dj4, dj2, dj6):
        return S.Zero
    if not _check_dj_couple(dj4, dj5, dj3):
        return S.Zero
    j123 = (dj1 + dj2 + dj3) // 2
    j156 = (dj1 + dj5 + dj6) // 2
    j426 = (dj4 + dj2 + dj6) // 2
    j453 = (dj4 + dj5 + dj3) // 2
    jpm123 = (dj1 + dj2 - dj3) // 2
    jpm132 = (dj1 + dj3 - dj2) // 2
    jpm231 = (dj2 + dj3 - dj1) // 2
    jpm156 = (dj1 + dj5 - dj6) // 2
    jpm426 = (dj4 + dj2 - dj6) // 2
    jpm453 = (dj4 + dj5 - dj3) // 2
    sumres = 0
    imin = max(j123, j453, j426, j156)
    imax = max(jpm123 + j453, jpm132 + j426, jpm231 + j156)
    for ii in range(int(imin), int(imax) + 1):
        ti = comb(ii + 1, j123 + 1) * \
           comb(jpm123, ii - j453) * \
           comb(jpm132, ii - j426) * \
           comb(jpm231, ii - j156)
        sumres = ti - sumres
    if sumres == 0:
        return S.Zero
    if imax % 2 == 1:
        sumres = -sumres
    res_num = comb(j123 + 1, dj1 + 1) * comb(dj1, jpm123)
    res_den = comb(j156 + 1, dj1 + 1) * comb(dj1, jpm156) * \
        comb(j426 + 1, dj4 + 1) * comb(dj4, jpm426) * \
        comb(j453 + 1, dj4 + 1) * comb(dj4, jpm453)
    res_den = res_den * (dj4 + 1) ** 2
    ressqrt = sqrt(Integer(res_num) / Integer(res_den))
    result = ressqrt * sumres
    if prec:
        result = result.evalf(prec)
    return result


def racah(aa, bb, cc, dd, ee, ff, prec=None):
    r"""
    Calculate the Racah symbol `W(a,b,c,d;e,f)`.

    Parameters
    ==========

    a, ..., f :
        Integer or half integer.
    prec :
        Precision, default: ``None``.

    Returns
    =======

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    Examples
    ========

    >>> from sympy.physics.wigner import racah
    >>> racah(3,3,3,3,3,3)
    -1/14

    Notes
    =====

    The Racah symbol is related to the Wigner 6j symbol:

    .. math::

       \operatorname{Wigner6j}(j_1,j_2,j_3,j_4,j_5,j_6)
       =(-1)^{j_1+j_2+j_4+j_5} W(j_1,j_2,j_5,j_4,j_3,j_6)

    Please see the 6j symbol for its much richer symmetries and for
    additional properties.

    Algorithm
    =========

    This function uses the algorithm of [Wei99]_ to calculate the
    value of the 6j symbol exactly. Note that the formula contains
    alternating sums over large binomials and is therefore unsuitable
    for finite precision arithmetic and only useful for a computer
    algebra system [Rasch03]_.

    Authors
    =======

    - Jens Rasch (2009-03-24): initial version
    """
    res = wigner_6j(aa, bb, ee, dd, cc, ff)
    if res == 0:
        return S.Zero
    phase = (-1) ** int(aa + bb + cc + dd)
    res = phase * res
    if prec:
        res = res.evalf(prec)
    return res


def wigner_9j(j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9, prec=None):
    r"""
    Calculate the Wigner 9j symbol
    `\operatorname{Wigner9j}(j_1,j_2,j_3,j_4,j_5,j_6,j_7,j_8,j_9)`.

    Parameters
    ==========

    j_1, ..., j_9 :
        Integer or half integer.
    prec :
        precision, default: ``None``.

    Returns
    =======

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    Examples
    ========

    >>> from sympy.physics.wigner import wigner_9j
    >>> wigner_9j(1,1,1, 1,1,1, 1,1,0, prec=64)
    0.05555555555555555555555555555555555555555555555555555555555555555

    >>> wigner_9j(1/2,1/2,0, 1/2,3/2,1, 0,1,1, prec=64)
    0.1666666666666666666666666666666666666666666666666666666666666667

    It is an error to have arguments that are not integer or half
    integer values::

        >>> wigner_9j(0.51,0.5,0.5, 0.5,0.5,0.5, 0.5,0.5,0.5)
        Traceback (most recent call last):
        ...
        ValueError: expecting integer or half-integer, got 0.51

    Algorithm
    =========

    This function uses the algorithm of [Wei99]_ to calculate the
    value of the 3j symbol exactly. Note that the formula contains
    alternating sums over large binomials and is therefore unsuitable
    for finite precision arithmetic and only useful for a computer
    algebra system [Rasch03]_.
    """
    dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9 = \
        map(_doubled_int, [j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9])
    if not _check_dj_couple(dj1, dj2, dj3):
        return S.Zero
    if not _check_dj_couple(dj4, dj5, dj6):
        return S.Zero
    if not _check_dj_couple(dj7, dj8, dj9):
        return S.Zero
    if not _check_dj_couple(dj1, dj4, dj7):
        return S.Zero
    if not _check_dj_couple(dj2, dj5, dj8):
        return S.Zero
    if not _check_dj_couple(dj3, dj6, dj9):
        return S.Zero
    j123 = (dj1 + dj2 + dj3) // 2
    j456 = (dj4 + dj5 + dj6) // 2
    j789 = (dj7 + dj8 + dj9) // 2
    j147 = (dj1 + dj4 + dj7) // 2
    j258 = (dj2 + dj5 + dj8) // 2
    j369 = (dj3 + dj6 + dj9) // 2
    pm123 = (dj1 + dj2 - dj3) // 2
    pm132 = (dj1 + dj3 - dj2) // 2
    pm231 = (dj2 + dj3 - dj1) // 2
    pm456 = (dj4 + dj5 - dj6) // 2
    pm465 = (dj4 + dj6 - dj5) // 2
    pm564 = (dj5 + dj6 - dj4) // 2
    pm789 = (dj7 + dj8 - dj9) // 2
    pm798 = (dj7 + dj9 - dj8) // 2
    pm897 = (dj8 + dj9 - dj7) // 2
    sumres = 0
    dtmin = max(abs(dj2 - dj6), abs(dj4 - dj8), abs(dj1 - dj9))
    dtmax = min(dj2 + dj6, dj4 + dj8, dj1 + dj9)
    for dt in range(dtmin, dtmax + 1, 2):
        j19t = (dj1 + dj9 + dt) // 2
        j26t = (dj2 + dj6 + dt) // 2
        j48t = (dj4 + dj8 + dt) // 2
        abc = dt + 1
        xmin = max(j123, j369, j26t, j19t)
        xmax = min(pm123 + j369, pm132 + j26t, pm231 + j19t)
        suma = 0
        for x in range(xmin, xmax + 1):
            ax = comb(x + 1, j19t + 1) * \
                comb((dj1 + dj9 - dt) // 2, x - j26t) * \
                comb((dj1 + dt - dj9) // 2, x - j369) * \
                comb((dt + dj9 - dj1) // 2, x - j123)
            suma = ax - suma
        ymin = max(j456, j258, j26t, j48t)
        ymax = min(pm456 + j26t, pm465 + j258, pm564 + j48t)
        sumb = 0
        for y in range(ymin, ymax + 1):
            by = comb(y + 1, j26t + 1) * \
                comb((dj2 + dj6 - dt) // 2, y - j48t) * \
                comb((dt + dj6 - dj2) // 2, y - j258) * \
                comb((dj2 + dt - dj6) // 2, y - j456)
            sumb = by - sumb
        zmin = max(j789, j147, j48t, j19t)
        zmax = min(pm789 + j19t, pm798 + j48t, pm897 + j147)
        sumc = 0
        for z in range(zmin, zmax + 1):
            cz = comb(z + 1, j48t + 1) * \
                comb((dj4 + dj8 - dt) // 2, z - j19t) * \
                comb((dt + dj8 - dj4) // 2, z - j147) * \
                comb((dj4 + dt - dj8) // 2, z - j789)
            sumc = cz - sumc
        abc = abc * suma * sumb * sumc
        if (xmax + ymax + zmax) % 2 == 1:
            abc = -abc
        sumres += abc
    if sumres == 0:
        return S.Zero
    if dtmax % 2 == 1:
        sumres = -sumres

    res_den = (dj3 + 1) * (dj6 + 1) * (dj9 + 1) * \
        (dj7 + 1) * (dj8 + 1) * (dj9 + 1)
    res_den *= comb(j123 + 1, dj3 + 1) * comb(dj3, pm231) * \
        comb(j456 + 1, dj6 + 1) * comb(dj6, pm564) * \
        comb(j789 + 1, dj9 + 1) * comb(dj9, pm897) * \
        comb(j147 + 1, dj7 + 1) * comb(dj7, (dj4 + dj7 - dj1) // 2) * \
        comb(j258 + 1, dj8 + 1) * comb(dj8, (dj5 + dj8 - dj2) // 2) * \
        comb(j369 + 1, dj9 + 1) * comb(dj9, (dj6 + dj9 - dj3) // 2)

    res = sumres * sqrt(1/Integer(res_den))
    if prec:
        res = res.evalf(prec)
    return res


def gaunt(l_1, l_2, l_3, m_1, m_2, m_3, prec=None):
    r"""
    Calculate the Gaunt coefficient.

    Explanation
    ===========

    The Gaunt coefficient is defined as the integral over three
    spherical harmonics:

    .. math::

        \begin{aligned}
        \operatorname{Gaunt}(l_1,l_2,l_3,m_1,m_2,m_3)
        &=\int Y_{l_1,m_1}(\Omega)
         Y_{l_2,m_2}(\Omega) Y_{l_3,m_3}(\Omega) \,d\Omega \\
        &=\sqrt{\frac{(2l_1+1)(2l_2+1)(2l_3+1)}{4\pi}}
         \operatorname{Wigner3j}(l_1,l_2,l_3,0,0,0)
         \operatorname{Wigner3j}(l_1,l_2,l_3,m_1,m_2,m_3)
        \end{aligned}

    Parameters
    ==========

    l_1, l_2, l_3, m_1, m_2, m_3 :
        Integer.
    prec:
        precision, default: ``None``.

    Returns
    =======

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    Examples
    ========

    >>> from sympy.physics.wigner import gaunt
    >>> gaunt(1,0,1,1,0,-1)
    -1/(2*sqrt(pi))
    >>> gaunt(1000,1000,1200,9,3,-12).n(64)
    0.006895004219221134484332976156744208248842039317638217822322799675

    It is an error to use non-integer values for `l` and `m`::

        >>> gaunt(1.2,0,1.2,0,0,0)
        Traceback (most recent call last):
        ...
        ValueError: 1.2 is not an integer
        >>> gaunt(1,0,1,1.1,0,-1.1)
        Traceback (most recent call last):
        ...
        ValueError: 1.1 is not an integer

    Notes
    =====

    The Gaunt coefficient obeys the following symmetry rules:

    - invariant under any permutation of the columns

      .. math::
        \begin{aligned}
          Y(l_1,l_2,l_3,m_1,m_2,m_3)
          &=Y(l_3,l_1,l_2,m_3,m_1,m_2) \\
          &=Y(l_2,l_3,l_1,m_2,m_3,m_1) \\
          &=Y(l_3,l_2,l_1,m_3,m_2,m_1) \\
          &=Y(l_1,l_3,l_2,m_1,m_3,m_2) \\
          &=Y(l_2,l_1,l_3,m_2,m_1,m_3)
        \end{aligned}

    - invariant under space inflection, i.e.

      .. math::
          Y(l_1,l_2,l_3,m_1,m_2,m_3)
          =Y(l_1,l_2,l_3,-m_1,-m_2,-m_3)

    - symmetric with respect to the 72 Regge symmetries as inherited
      for the `3j` symbols [Regge58]_

    - zero for `l_1`, `l_2`, `l_3` not fulfilling triangle relation

    - zero for violating any one of the conditions: `l_1 \ge |m_1|`,
      `l_2 \ge |m_2|`, `l_3 \ge |m_3|`

    - non-zero only for an even sum of the `l_i`, i.e.
      `L = l_1 + l_2 + l_3 = 2n` for `n` in `\mathbb{N}`

    Algorithms
    ==========

    This function uses the algorithm of [Wei99]_ and [Varshalovich88]_ to
    calculate the value of the Gaunt coefficient exactly. Note that
    the formula contains alternating sums over large binomials and is
    therefore unsuitable for finite precision arithmetic and only
    useful for a computer algebra system [Rasch03]_.

    Authors
    =======

    Jens Rasch (2009-03-24): initial version for Sage.
    """
    l1, l2, l3, m1, m2, m3 = map(_as_int, (l_1, l_2, l_3, m_1, m_2, m_3))

    if not _check_dj_couple(l1, l2, l3):
        return S.Zero
    if (m1 + m2 + m3) != 0:
        return S.Zero
    if (abs(m1) > l1) or (abs(m2) > l2) or (abs(m3) > l3):
        return S.Zero

    sumj = l1 + l2 + l3
    g = sumj // 2
    jm1 = l2 + l3 - l1
    jm2 = l1 + l3 - l2
    jm3 = l1 + l2 - l3
    j1mm1 = l1 - m1
    j2mm2 = l2 - m2
    j3mm3 = l3 - m3
    j1pm1 = l1 + m1
    imin = max(0, j1pm1 - jm2, j2mm2 - jm1)
    imax = min(jm3, j1pm1, j2mm2)

    sumres = 0
    for ii in range(int(imin), int(imax) + 1):
        ti = comb(jm3, ii) * \
           comb(jm2, j1pm1 - ii) * \
           comb(jm1, j2mm2 - ii)
        sumres = ti - sumres
    if sumres == 0:
        return S.Zero
    if ((g + j3mm3 + imax) % 2) == 1:
        sumres = -sumres

    sumres *= comb(g, g - l3) * comb(l3, g - l1)
    res_num = (2 * l3 + 1) * comb(2 * l3, jm1)
    res_den = comb(2 * l1, j1mm1) * \
        comb(2 * l2, j2mm2) * comb(2 * l3, j3mm3) * \
        comb(sumj + 1, 2 * l1 + 1) * comb(sumj + 1, 2 * l2 + 1)
    ressqrt = sqrt(Integer(res_num) / res_den)

    res = (ressqrt * sumres) / sqrt(4 * pi)
    if prec is not None:
        res = res.n(prec)
    return res


def real_gaunt(l_1, l_2, l_3, mu_1, mu_2, mu_3, prec=None):
    r"""
    Calculate the real Gaunt coefficient.

    Explanation
    ===========

    The real Gaunt coefficient is defined as the integral over three
    real spherical harmonics:

    .. math::
        \begin{aligned}
        \operatorname{RealGaunt}(l_1,l_2,l_3,\mu_1,\mu_2,\mu_3)
        &=\int Z^{\mu_1}_{l_1}(\Omega)
         Z^{\mu_2}_{l_2}(\Omega) Z^{\mu_3}_{l_3}(\Omega) \,d\Omega \\
        \end{aligned}

    Alternatively, it can be defined in terms of the standard Gaunt
    coefficient by relating the real spherical harmonics to the standard
    spherical harmonics via a unitary transformation `U`, i.e.
    `Z^{\mu}_{l}(\Omega)=\sum_{m'}U^{\mu}_{m'}Y^{m'}_{l}(\Omega)` [Homeier96]_.
    The real Gaunt coefficient is then defined as

    .. math::
        \begin{aligned}
        \operatorname{RealGaunt}(l_1,l_2,l_3,\mu_1,\mu_2,\mu_3)
        &=\int Z^{\mu_1}_{l_1}(\Omega)
         Z^{\mu_2}_{l_2}(\Omega) Z^{\mu_3}_{l_3}(\Omega) \,d\Omega \\
        &=\sum_{m'_1 m'_2 m'_3} U^{\mu_1}_{m'_1}U^{\mu_2}_{m'_2}U^{\mu_3}_{m'_3}
         \operatorname{Gaunt}(l_1,l_2,l_3,m'_1,m'_2,m'_3)
        \end{aligned}

    The unitary matrix `U` has components

    .. math::
        \begin{aligned}
        U^\mu_{m} = \delta_{|\mu||m|}*(\delta_{m0}\delta_{\mu 0} + \frac{1}{\sqrt{2}}\big[\Theta(\mu)\big(\delta_{m\mu}+(-1)^{m}\delta_{m-\mu}\big)
        +i \Theta(-\mu)\big((-1)^{m}\delta_{m\mu}-\delta_{m-\mu}\big)\big])
        \end{aligned}


    where `\delta_{ij}` is the Kronecker delta symbol and `\Theta` is a step
    function defined as

    .. math::
        \begin{aligned}
        \Theta(x) = \begin{cases} 1 \,\text{for}\, x > 0 \\ 0 \,\text{for}\, x \leq 0 \end{cases}
        \end{aligned}

    Parameters
    ==========

    l_1, l_2, l_3, mu_1, mu_2, mu_3 :
        Integer degree and order.

    prec:
        precision, default: ``None``.

    Returns
    =======

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    Examples
    ========
    >>> from sympy.physics.wigner import real_gaunt
    >>> real_gaunt(1,1,2,-1,1,-2)
    sqrt(15)/(10*sqrt(pi))
    >>> real_gaunt(10,10,20,-9,-9,0,prec=64)
    -0.00002480019791932209313156167176797577821140084216297395518482071448

    It is an error to use non-integer values for `l` and `\mu`::
        >>> real_gaunt(2.8,0.5,1.3,0,0,0)
        Traceback (most recent call last):
        ...
        ValueError: 2.8 is not an integer

        >>> real_gaunt(2,2,4,0.7,1,-3.4)
        Traceback (most recent call last):
        ...
        ValueError: 0.7 is not an integer

    Notes
    =====

    The real Gaunt coefficient inherits from the standard Gaunt coefficient,
    the invariance under any permutation of the pairs `(l_i, \mu_i)` and the
    requirement that the sum of the `l_i` be even to yield a non-zero value.
    It also obeys the following symmetry rules:

    - zero for `l_1`, `l_2`, `l_3` not fulfilling the condition
      `l_1 \in \{l_{\text{max}}, l_{\text{max}}-2, \ldots, l_{\text{min}}\}`,
      where `l_{\text{max}} = l_2+l_3`,

      .. math::
          \begin{aligned}
          l_{\text{min}} = \begin{cases} \kappa(l_2, l_3, \mu_2, \mu_3) & \text{if}\,
          \kappa(l_2, l_3, \mu_2, \mu_3) + l_{\text{max}}\, \text{is even} \\
          \kappa(l_2, l_3, \mu_2, \mu_3)+1 & \text{if}\, \kappa(l_2, l_3, \mu_2, \mu_3) +
          l_{\text{max}}\, \text{is odd}\end{cases}
          \end{aligned}

      and `\kappa(l_2, l_3, \mu_2, \mu_3) = \max{\big(|l_2-l_3|, \min{\big(|\mu_2+\mu_3|,
      |\mu_2-\mu_3|\big)}\big)}`

    - zero for an odd number of negative `\mu_i`

    Algorithms
    ==========

    This function uses the algorithms of [Homeier96]_ and [Wei99]_ to
    calculate the value of the real Gaunt coefficient exactly. Note that
    the formula used in [Wei99]_ contains alternating sums over large
    binomials and is therefore unsuitable for finite precision arithmetic
    and only useful for a computer algebra system [Wei99]_. However, this
    function can in principle use any algorithm that computes the Gaunt
    coefficient, so it is suitable for finite precision arithmetic in so far
    as the algorithm which computes the Gaunt coefficient is.
    """
    l_1, l_2, l_3, mu_1, mu_2, mu_3 = \
         map(_as_int, [l_1, l_2, l_3, mu_1, mu_2, mu_3])

    # check for quick exits
    if sum(1 for i in (mu_1, mu_2, mu_3) if i < 0) % 2:
        return S.Zero  # odd number of negative m
    if (l_1 + l_2 + l_3) % 2:
        return S.Zero  # sum of l is odd
    lmax = l_2 + l_3
    lmin = max(abs(l_2 - l_3), min(abs(mu_2 + mu_3), abs(mu_2 - mu_3)))
    if (lmin + lmax) % 2:
        lmin += 1
    if lmin not in range(lmax, lmin - 2, -2):
        return S.Zero

    kron_del = lambda i, j: 1 if i == j else 0
    s = lambda e: -1 if e % 2 else 1  #  (-1)**e to give +/-1, avoiding float when e<0

    t = lambda x: 1 if x > 0 else 0
    A = lambda mu, m: t(-mu) * (s(m) * kron_del(m, mu) - kron_del(m, -mu))
    B = lambda mu, m: t(mu) * (kron_del(m, mu) + s(m) * kron_del(m, -mu))
    U = lambda mu, m: kron_del(abs(mu), abs(m)) * (kron_del(mu, 0) * kron_del(m, 0) + (B(mu, m) + I * A(mu, m))/sqrt(2))

    ugnt = 0
    for m1 in range(-l_1, l_1+1):
        U1 = U(mu_1, m1)
        for m2 in range(-l_2, l_2+1):
            U2 = U(mu_2, m2)
            U3 = U(mu_3,-m1-m2)
            ugnt = ugnt + re(U1*U2*U3)*gaunt(l_1, l_2, l_3, m1, m2, -m1 - m2, prec=prec)

    return ugnt


class Wigner3j(Function):

    def doit(self, **hints):
        if all(obj.is_number for obj in self.args):
            return wigner_3j(*self.args)
        else:
            return self

def dot_rot_grad_Ynm(j, p, l, m, theta, phi):
    r"""
    Returns dot product of rotational gradients of spherical harmonics.

    Explanation
    ===========

    This function returns the right hand side of the following expression:

    .. math ::
        \vec{R}Y{_j^{p}} \cdot \vec{R}Y{_l^{m}} = (-1)^{m+p}
        \sum\limits_{k=|l-j|}^{l+j}Y{_k^{m+p}}  * \alpha_{l,m,j,p,k} *
        \frac{1}{2} (k^2-j^2-l^2+k-j-l)


    Arguments
    =========

    j, p, l, m .... indices in spherical harmonics (expressions or integers)
    theta, phi .... angle arguments in spherical harmonics

    Example
    =======

    >>> from sympy import symbols
    >>> from sympy.physics.wigner import dot_rot_grad_Ynm
    >>> theta, phi = symbols("theta phi")
    >>> dot_rot_grad_Ynm(3, 2, 2, 0, theta, phi).doit()
    3*sqrt(55)*Ynm(5, 2, theta, phi)/(11*sqrt(pi))

    """
    j = sympify(j)
    p = sympify(p)
    l = sympify(l)
    m = sympify(m)
    theta = sympify(theta)
    phi = sympify(phi)
    k = Dummy("k")

    def alpha(l,m,j,p,k):
        return sqrt((2*l+1)*(2*j+1)*(2*k+1)/(4*pi)) * \
                Wigner3j(j, l, k, S.Zero, S.Zero, S.Zero) * \
                Wigner3j(j, l, k, p, m, -m-p)

    return (S.NegativeOne)**(m+p) * Sum(Ynm(k, m+p, theta, phi) * alpha(l,m,j,p,k) / 2 \
        *(k**2-j**2-l**2+k-j-l), (k, abs(l-j), l+j))


def wigner_d_small(J, beta):
    """Return the small Wigner d matrix for angular momentum J.

    Explanation
    ===========

    J : An integer, half-integer, or SymPy symbol for the total angular
        momentum of the angular momentum space being rotated.
    beta : A real number representing the Euler angle of rotation about
        the so-called line of nodes. See [Edmonds74]_.

    Returns
    =======

    A matrix representing the corresponding Euler angle rotation( in the basis
    of eigenvectors of `J_z`).

    .. math ::
        \\mathcal{d}_{\\beta} = \\exp\\big( \\frac{i\\beta}{\\hbar} J_y\\big)

    such that

    .. math ::
        d^{(J)}_{m',m}(\\beta) = \\mathtt{wigner\\_d\\_small(J,beta)[J-mprime,J-m]}

    The components are calculated using the general form [Edmonds74]_,
    equation 4.1.15.

    Examples
    ========

    >>> from sympy import Integer, symbols, pi, pprint
    >>> from sympy.physics.wigner import wigner_d_small
    >>> half = 1/Integer(2)
    >>> beta = symbols("beta", real=True)
    >>> pprint(wigner_d_small(half, beta), use_unicode=True)
    ⎡   ⎛β⎞      ⎛β⎞⎤
    ⎢cos⎜─⎟   sin⎜─⎟⎥
    ⎢   ⎝2⎠      ⎝2⎠⎥
    ⎢               ⎥
    ⎢    ⎛β⎞     ⎛β⎞⎥
    ⎢-sin⎜─⎟  cos⎜─⎟⎥
    ⎣    ⎝2⎠     ⎝2⎠⎦

    >>> pprint(wigner_d_small(2*half, beta), use_unicode=True)
    ⎡        2⎛β⎞              ⎛β⎞    ⎛β⎞           2⎛β⎞     ⎤
    ⎢     cos ⎜─⎟        √2⋅sin⎜─⎟⋅cos⎜─⎟        sin ⎜─⎟     ⎥
    ⎢         ⎝2⎠              ⎝2⎠    ⎝2⎠            ⎝2⎠     ⎥
    ⎢                                                        ⎥
    ⎢       ⎛β⎞    ⎛β⎞       2⎛β⎞      2⎛β⎞        ⎛β⎞    ⎛β⎞⎥
    ⎢-√2⋅sin⎜─⎟⋅cos⎜─⎟  - sin ⎜─⎟ + cos ⎜─⎟  √2⋅sin⎜─⎟⋅cos⎜─⎟⎥
    ⎢       ⎝2⎠    ⎝2⎠        ⎝2⎠       ⎝2⎠        ⎝2⎠    ⎝2⎠⎥
    ⎢                                                        ⎥
    ⎢        2⎛β⎞               ⎛β⎞    ⎛β⎞          2⎛β⎞     ⎥
    ⎢     sin ⎜─⎟        -√2⋅sin⎜─⎟⋅cos⎜─⎟       cos ⎜─⎟     ⎥
    ⎣         ⎝2⎠               ⎝2⎠    ⎝2⎠           ⎝2⎠     ⎦

    From table 4 in [Edmonds74]_

    >>> pprint(wigner_d_small(half, beta).subs({beta:pi/2}), use_unicode=True)
    ⎡ √2   √2⎤
    ⎢ ──   ──⎥
    ⎢ 2    2 ⎥
    ⎢        ⎥
    ⎢-√2   √2⎥
    ⎢────  ──⎥
    ⎣ 2    2 ⎦

    >>> pprint(wigner_d_small(2*half, beta).subs({beta:pi/2}),
    ... use_unicode=True)
    ⎡       √2      ⎤
    ⎢1/2    ──   1/2⎥
    ⎢       2       ⎥
    ⎢               ⎥
    ⎢-√2         √2 ⎥
    ⎢────   0    ── ⎥
    ⎢ 2          2  ⎥
    ⎢               ⎥
    ⎢      -√2      ⎥
    ⎢1/2   ────  1/2⎥
    ⎣       2       ⎦

    >>> pprint(wigner_d_small(3*half, beta).subs({beta:pi/2}),
    ... use_unicode=True)
    ⎡ √2    √6    √6   √2⎤
    ⎢ ──    ──    ──   ──⎥
    ⎢ 4     4     4    4 ⎥
    ⎢                    ⎥
    ⎢-√6   -√2    √2   √6⎥
    ⎢────  ────   ──   ──⎥
    ⎢ 4     4     4    4 ⎥
    ⎢                    ⎥
    ⎢ √6   -√2   -√2   √6⎥
    ⎢ ──   ────  ────  ──⎥
    ⎢ 4     4     4    4 ⎥
    ⎢                    ⎥
    ⎢-√2    √6   -√6   √2⎥
    ⎢────   ──   ────  ──⎥
    ⎣ 4     4     4    4 ⎦

    >>> pprint(wigner_d_small(4*half, beta).subs({beta:pi/2}),
    ... use_unicode=True)
    ⎡             √6            ⎤
    ⎢1/4   1/2    ──   1/2   1/4⎥
    ⎢             4             ⎥
    ⎢                           ⎥
    ⎢-1/2  -1/2   0    1/2   1/2⎥
    ⎢                           ⎥
    ⎢ √6                     √6 ⎥
    ⎢ ──    0    -1/2   0    ── ⎥
    ⎢ 4                      4  ⎥
    ⎢                           ⎥
    ⎢-1/2  1/2    0    -1/2  1/2⎥
    ⎢                           ⎥
    ⎢             √6            ⎥
    ⎢1/4   -1/2   ──   -1/2  1/4⎥
    ⎣             4             ⎦

    """
    M = [J-i for i in range(2*J+1)]
    d = zeros(2*J+1)

    # Mi corresponds to Edmonds' $m'$, and Mj to $m$.
    for i, Mi in enumerate(M):
        for j, Mj in enumerate(M):

            # We get the maximum and minimum value of sigma.
            sigmamax = min([J-Mi, J-Mj])
            sigmamin = max([0, -Mi-Mj])

            dij = sqrt(binomial(2*J, J+Mj) /
                       binomial(2*J, J+Mi))
            terms = [(-1)**(J-Mi-s) *
                     binomial(J+Mj, J-Mi-s) *
                     binomial(J-Mj, s) *
                     cos(beta/2)**(2*s+Mi+Mj) *
                     sin(beta/2)**(2*J-2*s-Mj-Mi)
                     for s in range(sigmamin, sigmamax+1)]

            d[i, j] = dij*Add(*terms)

    return ImmutableMatrix(d)


def wigner_d(J, alpha, beta, gamma):
    """Return the Wigner D matrix for angular momentum J.

    Explanation
    ===========

    J :
        An integer, half-integer, or SymPy symbol for the total angular
        momentum of the angular momentum space being rotated.
    alpha, beta, gamma - Real numbers representing the Euler.
        Angles of rotation about the so-called figure axis, line of nodes,
        and vertical. See [Edmonds74]_, however note that the symbols alpha
        and gamma are swapped in this implementation.

    Returns
    =======

    A matrix representing the corresponding Euler angle rotation (in the basis
    of eigenvectors of `J_z`).

    .. math ::
        \\mathcal{D}_{\\alpha \\beta \\gamma} =
        \\exp\\big( \\frac{i\\alpha}{\\hbar} J_z\\big)
        \\exp\\big( \\frac{i\\beta}{\\hbar} J_y\\big)
        \\exp\\big( \\frac{i\\gamma}{\\hbar} J_z\\big)

    such that

    .. math ::
        \\mathcal{D}^{(J)}_{m',m}(\\alpha, \\beta, \\gamma) =
        \\mathtt{wigner_d(J, alpha, beta, gamma)[J-mprime,J-m]}

    The components are calculated using the general form [Edmonds74]_,
    equation 4.1.12, however note that the angles alpha and gamma are swapped
    in this implementation.

    Examples
    ========

    The simplest possible example:

    >>> from sympy.physics.wigner import wigner_d
    >>> from sympy import Integer, symbols, pprint
    >>> half = 1/Integer(2)
    >>> alpha, beta, gamma = symbols("alpha, beta, gamma", real=True)
    >>> pprint(wigner_d(half, alpha, beta, gamma), use_unicode=True)
    ⎡  ⅈ⋅α  ⅈ⋅γ             ⅈ⋅α  -ⅈ⋅γ         ⎤
    ⎢  ───  ───             ───  ─────        ⎥
    ⎢   2    2     ⎛β⎞       2     2      ⎛β⎞ ⎥
    ⎢ ℯ   ⋅ℯ   ⋅cos⎜─⎟     ℯ   ⋅ℯ     ⋅sin⎜─⎟ ⎥
    ⎢              ⎝2⎠                    ⎝2⎠ ⎥
    ⎢                                         ⎥
    ⎢  -ⅈ⋅α   ⅈ⋅γ          -ⅈ⋅α   -ⅈ⋅γ        ⎥
    ⎢  ─────  ───          ─────  ─────       ⎥
    ⎢    2     2     ⎛β⎞     2      2      ⎛β⎞⎥
    ⎢-ℯ     ⋅ℯ   ⋅sin⎜─⎟  ℯ     ⋅ℯ     ⋅cos⎜─⎟⎥
    ⎣                ⎝2⎠                   ⎝2⎠⎦

    """
    d = wigner_d_small(J, beta)
    M = [J-i for i in range(2*J+1)]
    # Mi corresponds to Edmonds' $m'$, and Mj to $m$.
    D = [[exp(I*Mi*alpha)*d[i, j]*exp(I*Mj*gamma)
          for j, Mj in enumerate(M)] for i, Mi in enumerate(M)]
    return ImmutableMatrix(D)
