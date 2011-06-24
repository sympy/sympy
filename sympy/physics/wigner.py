r"""
Wigner, Clebsch-Gordan, Racah, and Gaunt coefficients

Collection of functions for calculating Wigner 3j, 6j, 9j,
Clebsch-Gordan, Racah as well as Gaunt coefficients exactly, all
evaluating to a rational number times the square root of a rational
number [Rasch03]_.

Please see the description of the individual functions for further
details and examples.

REFERENCES:

.. [Rasch03] J. Rasch and A. C. H. Yu, 'Efficient Storage Scheme for
  Pre-calculated Wigner 3j, 6j and Gaunt Coefficients', SIAM
  J. Sci. Comput. Volume 25, Issue 4, pp. 1416-1428 (2003)

This code was taken from Sage with the permission of all authors:

http://groups.google.com/group/sage-devel/browse_thread/thread/33835976efbb3b7f

AUTHORS:

- Jens Rasch (2009-03-24): initial version for Sage

- Jens Rasch (2009-05-31): updated to sage-4.0
Copyright (C) 2008 Jens Rasch <jyr2000@gmail.com>
"""
from sympy import Integer, pi, sqrt
#from sage.rings.complex_number import ComplexNumber
#from sage.rings.finite_rings.integer_mod import Mod

# This list of precomputed factorials is needed to massively
# accelerate future calculations of the various coefficients
_Factlist=[1]

def _calc_factlist(nn):
    r"""
    Function calculates a list of precomputed factorials in order to
    massively accelerate future calculations of the various
    coefficients.

    INPUT:

    -  ``nn`` -  integer, highest factorial to be computed

    OUTPUT:

    list of integers -- the list of precomputed factorials

    EXAMPLES:

    Calculate list of factorials::

        sage: from sage.functions.wigner import _calc_factlist
        sage: _calc_factlist(10)
        [1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800]
    """
    if nn >= len(_Factlist):
        for ii in range(len(_Factlist), nn + 1):
            _Factlist.append(_Factlist[ii - 1] * ii)
    return _Factlist[:int(nn) + 1]


def wigner_3j(j_1, j_2, j_3, m_1, m_2, m_3, prec=None):
    r"""
    Calculate the Wigner 3j symbol `Wigner3j(j_1,j_2,j_3,m_1,m_2,m_3)`.

    INPUT:

    -  ``j_1``, ``j_2``, ``j_3``, ``m_1``, ``m_2``, ``m_3`` - integer or half integer

    -  ``prec`` - precision, default: ``None``. Providing a precision can
       drastically speed up the calculation.

    OUTPUT:

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    EXAMPLES::

        sage: wigner_3j(2, 6, 4, 0, 0, 0)
        sqrt(5/143)
        sage: wigner_3j(2, 6, 4, 0, 0, 1)
        0
        sage: wigner_3j(0.5, 0.5, 1, 0.5, -0.5, 0)
        sqrt(1/6)
        sage: wigner_3j(40, 100, 60, -10, 60, -50)
        95608/18702538494885*sqrt(21082735836735314343364163310/220491455010479533763)
        sage: wigner_3j(2500, 2500, 5000, 2488, 2400, -4888, prec=64)
        7.60424456883448589e-12

    It is an error to have arguments that are not integer or half
    integer values::

        sage: wigner_3j(2.1, 6, 4, 0, 0, 0)
        Traceback (most recent call last):
        ...
        ValueError: j values must be integer or half integer
        sage: wigner_3j(2, 6, 4, 1, 0, -1.1)
        Traceback (most recent call last):
        ...
        ValueError: m values must be integer or half integer

    NOTES:

    The Wigner 3j symbol obeys the following symmetry rules:

    - invariant under any permutation of the columns (with the
      exception of a sign change where `J:=j_1+j_2+j_3`):

      .. math::

         Wigner3j(j_1,j_2,j_3,m_1,m_2,m_3)
          =Wigner3j(j_3,j_1,j_2,m_3,m_1,m_2)
          =Wigner3j(j_2,j_3,j_1,m_2,m_3,m_1)
          =(-1)^J Wigner3j(j_3,j_2,j_1,m_3,m_2,m_1)
          =(-1)^J Wigner3j(j_1,j_3,j_2,m_1,m_3,m_2)
          =(-1)^J Wigner3j(j_2,j_1,j_3,m_2,m_1,m_3)

    - invariant under space inflection, i.e.

      .. math::

         Wigner3j(j_1,j_2,j_3,m_1,m_2,m_3)
         =(-1)^J Wigner3j(j_1,j_2,j_3,-m_1,-m_2,-m_3)

    - symmetric with respect to the 72 additional symmetries based on
      the work by [Regge58]_

    - zero for `j_1`, `j_2`, `j_3` not fulfilling triangle relation

    - zero for `m_1 + m_2 + m_3 \neq 0`

    - zero for violating any one of the conditions
      `j_1 \ge |m_1|`,  `j_2 \ge |m_2|`,  `j_3 \ge |m_3|`

    ALGORITHM:

    This function uses the algorithm of [Edmonds74]_ to calculate the
    value of the 3j symbol exactly. Note that the formula contains
    alternating sums over large factorials and is therefore unsuitable
    for finite precision arithmetic and only useful for a computer
    algebra system [Rasch03]_.

    REFERENCES:

    .. [Regge58] 'Symmetry Properties of Clebsch-Gordan Coefficients',
      T. Regge, Nuovo Cimento, Volume 10, pp. 544 (1958)

    .. [Edmonds74] 'Angular Momentum in Quantum Mechanics',
      A. R. Edmonds, Princeton University Press (1974)

    AUTHORS:

    - Jens Rasch (2009-03-24): initial version
    """
    if int(j_1 * 2) != j_1 * 2 or int(j_2 * 2) != j_2 * 2 or \
            int(j_3 * 2) != j_3 * 2:
        raise ValueError("j values must be integer or half integer")
    if int(m_1 * 2) != m_1 * 2 or int(m_2 * 2) != m_2 * 2 or \
            int(m_3 * 2) != m_3 * 2:
        raise ValueError("m values must be integer or half integer")
    if m_1 + m_2 + m_3 != 0:
        return 0
    prefid = Integer((-1) ** int(j_1 - j_2 - m_3))
    m_3 = -m_3
    a1 = j_1 + j_2 - j_3
    if a1 < 0:
        return 0
    a2 = j_1 - j_2 + j_3
    if a2 < 0:
        return 0
    a3 = -j_1 + j_2 + j_3
    if a3 < 0:
        return 0
    if (abs(m_1) > j_1) or (abs(m_2) > j_2) or (abs(m_3) > j_3):
        return 0

    maxfact = max(j_1 + j_2 + j_3 + 1, j_1 + abs(m_1), j_2 + abs(m_2), \
                  j_3 + abs(m_3))
    _calc_factlist(maxfact)

    argsqrt = Integer(_Factlist[int(j_1 + j_2 - j_3)] * \
                          _Factlist[int(j_1 - j_2 + j_3)] * \
                          _Factlist[int(-j_1 + j_2 + j_3)] * \
                          _Factlist[int(j_1 - m_1)] * \
                          _Factlist[int(j_1 + m_1)] * \
                          _Factlist[int(j_2 - m_2)] * \
                          _Factlist[int(j_2 + m_2)] * \
                          _Factlist[int(j_3 - m_3)] * \
                          _Factlist[int(j_3 + m_3)]) / \
                          _Factlist[int(j_1 + j_2 + j_3 + 1)]

    ressqrt = sqrt(argsqrt)
    if ressqrt.is_complex:
        ressqrt = ressqrt.as_real_imag()[0]

    imin = max(-j_3 + j_1 + m_2, -j_3 + j_2 - m_1, 0)
    imax = min(j_2 + m_2, j_1 - m_1, j_1 + j_2 - j_3)
    sumres = 0
    for ii in range(imin, imax + 1):
        den = _Factlist[ii] * \
            _Factlist[int(ii + j_3 - j_1 - m_2)] * \
            _Factlist[int(j_2 + m_2 - ii)] * \
            _Factlist[int(j_1 - ii - m_1)] * \
            _Factlist[int(ii + j_3 - j_2 + m_1)] * \
            _Factlist[int(j_1 + j_2 - j_3 - ii)]
        sumres = sumres + Integer((-1) ** ii) / den

    res = ressqrt * sumres * prefid
    return res


def clebsch_gordan(j_1, j_2, j_3, m_1, m_2, m_3, prec=None):
    r"""
    Calculates the Clebsch-Gordan coefficient
    `\langle j_1 m_1 \; j_2 m_2 | j_3 m_3 \rangle`.

    The reference for this function is [Edmonds74]_.

    INPUT:

    -  ``j_1``, ``j_2``, ``j_3``, ``m_1``, ``m_2``, ``m_3`` - integer or half integer

    -  ``prec`` - precision, default: ``None``. Providing a precision can
       drastically speed up the calculation.

    OUTPUT:

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    EXAMPLES::

        >>> from sympy import S
        >>> from sympy.physics.wigner import clebsch_gordan
        >>> clebsch_gordan(S(3)/2, S(1)/2, 2, S(3)/2, S(1)/2, 2)
        1
        >>> clebsch_gordan(S(3)/2, S(1)/2, 1, S(3)/2, -S(1)/2, 1)
        3**(1/2)/2
        >>> clebsch_gordan(S(3)/2, S(1)/2, 1, -S(1)/2, S(1)/2, 0)
        -2**(1/2)/2

    NOTES:

    The Clebsch-Gordan coefficient will be evaluated via its relation
    to Wigner 3j symbols:

    .. math::

        \langle j_1 m_1 \; j_2 m_2 | j_3 m_3 \rangle
        =(-1)^{j_1-j_2+m_3} \sqrt{2j_3+1} \;
        Wigner3j(j_1,j_2,j_3,m_1,m_2,-m_3)

    See also the documentation on Wigner 3j symbols which exhibit much
    higher symmetry relations than the Clebsch-Gordan coefficient.

    AUTHORS:

    - Jens Rasch (2009-03-24): initial version
    """
    res = (-1) ** int(j_1 - j_2 + m_3) * sqrt(2 * j_3 + 1) * \
        wigner_3j(j_1, j_2, j_3, m_1, m_2, -m_3, prec)
    return res


def _big_delta_coeff(aa, bb, cc, prec=None):
    r"""
    Calculates the Delta coefficient of the 3 angular momenta for
    Racah symbols. Also checks that the differences are of integer
    value.

    INPUT:

    -  ``aa`` - first angular momentum, integer or half integer

    -  ``bb`` - second angular momentum, integer or half integer

    -  ``cc`` - third angular momentum, integer or half integer

    -  ``prec`` - precision of the ``sqrt()`` calculation

    OUTPUT:

    double - Value of the Delta coefficient

    EXAMPLES::

        sage: from sage.functions.wigner import _big_delta_coeff
        sage: _big_delta_coeff(1,1,1)
        1/2*sqrt(1/6)
    """
    if int(aa + bb - cc) != (aa + bb - cc):
        raise ValueError("j values must be integer or half integer and fulfill the triangle relation")
    if int(aa + cc - bb) != (aa + cc - bb):
        raise ValueError("j values must be integer or half integer and fulfill the triangle relation")
    if int(bb + cc - aa) != (bb + cc - aa):
        raise ValueError("j values must be integer or half integer and fulfill the triangle relation")
    if (aa + bb - cc) < 0:
        return 0
    if (aa + cc - bb) < 0:
        return 0
    if (bb + cc - aa) < 0:
        return 0

    maxfact = max(aa + bb - cc, aa + cc - bb, bb + cc - aa, aa + bb + cc + 1)
    _calc_factlist(maxfact)

    argsqrt = Integer(_Factlist[int(aa + bb - cc)] * \
                          _Factlist[int(aa + cc - bb)] * \
                          _Factlist[int(bb + cc - aa)]) / \
                          Integer(_Factlist[int(aa + bb + cc + 1)])

    ressqrt = sqrt(argsqrt)
    if prec:
        ressqrt = ressqrt.evalf(prec).as_real_imag()[0]
    return ressqrt


def racah(aa, bb, cc, dd, ee, ff, prec=None):
    r"""
    Calculate the Racah symbol `W(a,b,c,d;e,f)`.

    INPUT:

    -  ``a``, ..., ``f`` - integer or half integer

    -  ``prec`` - precision, default: ``None``. Providing a precision can
       drastically speed up the calculation.

    OUTPUT:

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    EXAMPLES::

        sage: racah(3,3,3,3,3,3)
        -1/14

    NOTES:

    The Racah symbol is related to the Wigner 6j symbol:

    .. math::

       Wigner6j(j_1,j_2,j_3,j_4,j_5,j_6)
       =(-1)^{j_1+j_2+j_4+j_5} W(j_1,j_2,j_5,j_4,j_3,j_6)

    Please see the 6j symbol for its much richer symmetries and for
    additional properties.

    ALGORITHM:

    This function uses the algorithm of [Edmonds74]_ to calculate the
    value of the 6j symbol exactly. Note that the formula contains
    alternating sums over large factorials and is therefore unsuitable
    for finite precision arithmetic and only useful for a computer
    algebra system [Rasch03]_.

    AUTHORS:

    - Jens Rasch (2009-03-24): initial version
    """
    prefac = _big_delta_coeff(aa, bb, ee, prec) * \
        _big_delta_coeff(cc, dd, ee, prec) * \
        _big_delta_coeff(aa, cc, ff, prec) * \
        _big_delta_coeff(bb, dd, ff, prec)
    if prefac == 0:
        return 0
    imin = max(aa + bb + ee, cc + dd + ee, aa + cc + ff, bb + dd + ff)
    imax = min(aa + bb + cc + dd, aa + dd + ee + ff, bb + cc + ee + ff)

    maxfact = max(imax + 1, aa + bb + cc + dd, aa + dd + ee + ff, \
                      bb + cc + ee + ff)
    _calc_factlist(maxfact)

    sumres = 0
    for kk in range(imin, imax + 1):
        den = _Factlist[int(kk - aa - bb - ee)] * \
            _Factlist[int(kk - cc - dd - ee)] * \
            _Factlist[int(kk - aa - cc - ff)] * \
            _Factlist[int(kk - bb - dd - ff)] * \
            _Factlist[int(aa + bb + cc + dd - kk)] * \
            _Factlist[int(aa + dd + ee + ff - kk)] * \
            _Factlist[int(bb + cc + ee + ff - kk)]
        sumres = sumres + Integer((-1) ** kk * _Factlist[kk + 1]) / den

    res = prefac * sumres * (-1) ** int(aa + bb + cc + dd)
    return res


def wigner_6j(j_1, j_2, j_3, j_4, j_5, j_6, prec=None):
    r"""
    Calculate the Wigner 6j symbol `Wigner6j(j_1,j_2,j_3,j_4,j_5,j_6)`.

    INPUT:

    -  ``j_1``, ..., ``j_6`` - integer or half integer

    -  ``prec`` - precision, default: ``None``. Providing a precision can
       drastically speed up the calculation.

    OUTPUT:

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    EXAMPLES::

        sage: wigner_6j(3,3,3,3,3,3)
        -1/14
        sage: wigner_6j(5,5,5,5,5,5)
        1/52
        sage: wigner_6j(6,6,6,6,6,6)
        309/10868
        sage: wigner_6j(8,8,8,8,8,8)
        -12219/965770
        sage: wigner_6j(30,30,30,30,30,30)
        36082186869033479581/87954851694828981714124
        sage: wigner_6j(0.5,0.5,1,0.5,0.5,1)
        1/6
        sage: wigner_6j(200,200,200,200,200,200, prec=1000)*1.0
        0.000155903212413242

    It is an error to have arguments that are not integer or half
    integer values or do not fulfill the triangle relation::

        sage: wigner_6j(2.5,2.5,2.5,2.5,2.5,2.5)
        Traceback (most recent call last):
        ...
        ValueError: j values must be integer or half integer and fulfill the triangle relation
        sage: wigner_6j(0.5,0.5,1.1,0.5,0.5,1.1)
        Traceback (most recent call last):
        ...
        ValueError: j values must be integer or half integer and fulfill the triangle relation

    NOTES:

    The Wigner 6j symbol is related to the Racah symbol but exhibits
    more symmetries as detailed below.

    .. math::

       Wigner6j(j_1,j_2,j_3,j_4,j_5,j_6)
        =(-1)^{j_1+j_2+j_4+j_5} W(j_1,j_2,j_5,j_4,j_3,j_6)

    The Wigner 6j symbol obeys the following symmetry rules:

    - Wigner 6j symbols are left invariant under any permutation of
      the columns:

      .. math::

         Wigner6j(j_1,j_2,j_3,j_4,j_5,j_6)
          =Wigner6j(j_3,j_1,j_2,j_6,j_4,j_5)
          =Wigner6j(j_2,j_3,j_1,j_5,j_6,j_4)
          =Wigner6j(j_3,j_2,j_1,j_6,j_5,j_4)
          =Wigner6j(j_1,j_3,j_2,j_4,j_6,j_5)
          =Wigner6j(j_2,j_1,j_3,j_5,j_4,j_6)

    - They are invariant under the exchange of the upper and lower
      arguments in each of any two columns, i.e.

      .. math::

         Wigner6j(j_1,j_2,j_3,j_4,j_5,j_6)
          =Wigner6j(j_1,j_5,j_6,j_4,j_2,j_3)
          =Wigner6j(j_4,j_2,j_6,j_1,j_5,j_3)
          =Wigner6j(j_4,j_5,j_3,j_1,j_2,j_6)

    - additional 6 symmetries [Regge59]_ giving rise to 144 symmetries
      in total

    - only non-zero if any triple of `j`'s fulfill a triangle relation

    ALGORITHM:

    This function uses the algorithm of [Edmonds74]_ to calculate the
    value of the 6j symbol exactly. Note that the formula contains
    alternating sums over large factorials and is therefore unsuitable
    for finite precision arithmetic and only useful for a computer
    algebra system [Rasch03]_.

    REFERENCES:

    .. [Regge59] 'Symmetry Properties of Racah Coefficients',
      T. Regge, Nuovo Cimento, Volume 11, pp. 116 (1959)
    """
    res = (-1) ** int(j_1 + j_2 + j_4 + j_5) * \
        racah(j_1, j_2, j_5, j_4, j_3, j_6, prec)
    return res


def wigner_9j(j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9, prec=None):
    r"""
    Calculate the Wigner 9j symbol
    `Wigner9j(j_1,j_2,j_3,j_4,j_5,j_6,j_7,j_8,j_9)`.

    INPUT:

    -  ``j_1``, ..., ``j_9`` - integer or half integer

    -  ``prec`` - precision, default: ``None``. Providing a precision can
       drastically speed up the calculation.

    OUTPUT:

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    EXAMPLES:

    A couple of examples and test cases, note that for speed reasons a
    precision is given::

        sage: wigner_9j(1,1,1, 1,1,1, 1,1,0 ,prec=64) # ==1/18
        0.0555555555555555555
        sage: wigner_9j(1,1,1, 1,1,1, 1,1,1)
        0
        sage: wigner_9j(1,1,1, 1,1,1, 1,1,2 ,prec=64) # ==1/18
        0.0555555555555555556
        sage: wigner_9j(1,2,1, 2,2,2, 1,2,1 ,prec=64) # ==-1/150
        -0.00666666666666666667
        sage: wigner_9j(3,3,2, 2,2,2, 3,3,2 ,prec=64) # ==157/14700
        0.0106802721088435374
        sage: wigner_9j(3,3,2, 3,3,2, 3,3,2 ,prec=64) # ==3221*sqrt(70)/(246960*sqrt(105)) - 365/(3528*sqrt(70)*sqrt(105))
        0.00944247746651111739
        sage: wigner_9j(3,3,1, 3.5,3.5,2, 3.5,3.5,1 ,prec=64) # ==3221*sqrt(70)/(246960*sqrt(105)) - 365/(3528*sqrt(70)*sqrt(105))
        0.0110216678544351364
        sage: wigner_9j(100,80,50, 50,100,70, 60,50,100 ,prec=1000)*1.0
        1.05597798065761e-7
        sage: wigner_9j(30,30,10, 30.5,30.5,20, 30.5,30.5,10 ,prec=1000)*1.0 # ==(80944680186359968990/95103769817469)*sqrt(1/682288158959699477295)
        0.0000325841699408828
        sage: wigner_9j(64,62.5,114.5, 61.5,61,112.5, 113.5,110.5,60, prec=1000)*1.0
        -3.41407910055520e-39
        sage: wigner_9j(15,15,15, 15,3,15, 15,18,10, prec=1000)*1.0
        -0.0000778324615309539
        sage: wigner_9j(1.5,1,1.5, 1,1,1, 1.5,1,1.5)
        0

    It is an error to have arguments that are not integer or half
    integer values or do not fulfill the triangle relation::

        sage: wigner_9j(0.5,0.5,0.5, 0.5,0.5,0.5, 0.5,0.5,0.5,prec=64)
        Traceback (most recent call last):
        ...
        ValueError: j values must be integer or half integer and fulfill the triangle relation
        sage: wigner_9j(1,1,1, 0.5,1,1.5, 0.5,1,2.5,prec=64)
        Traceback (most recent call last):
        ...
        ValueError: j values must be integer or half integer and fulfill the triangle relation

    ALGORITHM:

    This function uses the algorithm of [Edmonds74]_ to calculate the
    value of the 3j symbol exactly. Note that the formula contains
    alternating sums over large factorials and is therefore unsuitable
    for finite precision arithmetic and only useful for a computer
    algebra system [Rasch03]_.
    """
    imin = 0
    imax = min(j_1 + j_9, j_2 + j_6, j_4 + j_8)

    sumres = 0
    for kk in range(imin, imax + 1):
        sumres = sumres + (2 * kk + 1) * \
            racah(j_1, j_2, j_9, j_6, j_3, kk, prec) * \
            racah(j_4, j_6, j_8, j_2, j_5, kk, prec) * \
            racah(j_1, j_4, j_9, j_8, j_7, kk, prec)
    return sumres


def gaunt(l_1, l_2, l_3, m_1, m_2, m_3, prec=None):
    r"""
    Calculate the Gaunt coefficient.

    The Gaunt coefficient is defined as the integral over three
    spherical harmonics:

    .. math::

        Y(j_1,j_2,j_3,m_1,m_2,m_3)
        =\int Y_{l_1,m_1}(\Omega)
         Y_{l_2,m_2}(\Omega) Y_{l_3,m_3}(\Omega) d\Omega
        =\sqrt{(2l_1+1)(2l_2+1)(2l_3+1)/(4\pi)}
         \; Y(j_1,j_2,j_3,0,0,0) \; Y(j_1,j_2,j_3,m_1,m_2,m_3)

    INPUT:

    -  ``l_1``, ``l_2``, ``l_3``, ``m_1``, ``m_2``, ``m_3`` - integer

    -  ``prec`` - precision, default: ``None``. Providing a precision can
       drastically speed up the calculation.

    OUTPUT:

    Rational number times the square root of a rational number
    (if ``prec=None``), or real number if a precision is given.

    EXAMPLES::

        sage: gaunt(1,0,1,1,0,-1)
        -1/2/sqrt(pi)
        sage: gaunt(1,0,1,1,0,0)
        0
        sage: gaunt(29,29,34,10,-5,-5)
        1821867940156/215552371055153321*sqrt(22134)/sqrt(pi)
        sage: gaunt(20,20,40,1,-1,0)
        28384503878959800/74029560764440771/sqrt(pi)
        sage: gaunt(12,15,5,2,3,-5)
        91/124062*sqrt(36890)/sqrt(pi)
        sage: gaunt(10,10,12,9,3,-12)
        -98/62031*sqrt(6279)/sqrt(pi)
        sage: gaunt(1000,1000,1200,9,3,-12).n(64)
        0.00689500421922113448

    It is an error to use non-integer values for `l` and `m`::

        sage: gaunt(1.2,0,1.2,0,0,0)
        Traceback (most recent call last):
        ...
        ValueError: l values must be integer
        sage: gaunt(1,0,1,1.1,0,-1.1)
        Traceback (most recent call last):
        ...
        ValueError: m values must be integer

    NOTES:

    The Gaunt coefficient obeys the following symmetry rules:

    - invariant under any permutation of the columns

      .. math::
          Y(j_1,j_2,j_3,m_1,m_2,m_3)
          =Y(j_3,j_1,j_2,m_3,m_1,m_2)
          =Y(j_2,j_3,j_1,m_2,m_3,m_1)
          =Y(j_3,j_2,j_1,m_3,m_2,m_1)
          =Y(j_1,j_3,j_2,m_1,m_3,m_2)
          =Y(j_2,j_1,j_3,m_2,m_1,m_3)

    - invariant under space inflection, i.e.

      .. math::
          Y(j_1,j_2,j_3,m_1,m_2,m_3)
          =Y(j_1,j_2,j_3,-m_1,-m_2,-m_3)

    - symmetric with respect to the 72 Regge symmetries as inherited
      for the `3j` symbols [Regge58]_

    - zero for `l_1`, `l_2`, `l_3` not fulfilling triangle relation

    - zero for violating any one of the conditions: `l_1 \ge |m_1|`,
      `l_2 \ge |m_2|`, `l_3 \ge |m_3|`

    - non-zero only for an even sum of the `l_i`, i.e.
      `J=l_1+l_2+l_3=2n` for `n` in `\Bold{N}`

    ALGORITHM:

    This function uses the algorithm of [Liberatodebrito82]_ to
    calculate the value of the Gaunt coefficient exactly. Note that
    the formula contains alternating sums over large factorials and is
    therefore unsuitable for finite precision arithmetic and only
    useful for a computer algebra system [Rasch03]_.

    REFERENCES:

    .. [Liberatodebrito82] 'FORTRAN program for the integral of three
      spherical harmonics', A. Liberato de Brito,
      Comput. Phys. Commun., Volume 25, pp. 81-85 (1982)

    AUTHORS:

    - Jens Rasch (2009-03-24): initial version for Sage
    """
    if int(l_1) != l_1 or int(l_2) != l_2 or int(l_3) != l_3:
        raise ValueError("l values must be integer")
    if int(m_1) != m_1 or int(m_2) != m_2 or int(m_3) != m_3:
        raise ValueError("m values must be integer")

    bigL = (l_1 + l_2 + l_3) // 2
    a1 = l_1 + l_2 - l_3
    if a1 < 0:
        return 0
    a2 = l_1 - l_2 + l_3
    if a2 < 0:
        return 0
    a3 = -l_1 + l_2 + l_3
    if a3 < 0:
        return 0
    if (2 * bigL) % 2 != 0:
        return 0
    if (m_1 + m_2 + m_3) != 0:
        return 0
    if (abs(m_1) > l_1) or (abs(m_2) > l_2) or (abs(m_3) > l_3):
        return 0

    imin = max(-l_3 + l_1 + m_2, -l_3 + l_2 - m_1, 0)
    imax = min(l_2 + m_2, l_1 - m_1, l_1 + l_2 - l_3)

    maxfact = max(l_1 + l_2 + l_3 + 1, imax + 1)
    _calc_factlist(maxfact)

    argsqrt = (2 * l_1 + 1) * (2 * l_2 + 1) * (2 * l_3 + 1) * \
        _Factlist[l_1 - m_1] * _Factlist[l_1 + m_1] * _Factlist[l_2 - m_2] * \
        _Factlist[l_2 + m_2] * _Factlist[l_3 - m_3] * _Factlist[l_3 + m_3] / \
        (4*pi)
    ressqrt = sqrt(argsqrt)

    prefac = Integer(_Factlist[bigL] * _Factlist[l_2 - l_1 + l_3] * \
                     _Factlist[l_1 - l_2 + l_3] * _Factlist[l_1 + l_2 - l_3])/ \
                     _Factlist[2 * bigL+1]/ \
                     (_Factlist[bigL - l_1] * _Factlist[bigL - l_2] * _Factlist[bigL - l_3])

    sumres = 0
    for ii in range(imin, imax + 1):
        den = _Factlist[ii] * _Factlist[ii + l_3 - l_1 - m_2] * \
            _Factlist[l_2 + m_2 - ii] * _Factlist[l_1 - ii - m_1] * \
            _Factlist[ii + l_3 - l_2 + m_1] * _Factlist[l_1 + l_2 - l_3 - ii]
        sumres = sumres + Integer((-1) ** ii) / den

    res = ressqrt * prefac * sumres * (-1) ** (bigL + l_3 + m_1 - m_2)
    if prec != None:
        res = res.n(prec)
    return res
