from sympy.core.compatibility import as_int
from sympy.core.numbers import Rational
from sympy.core.numbers import Integer

from fractions import gcd


def continued_fraction(num, den, delta):
    """
    Return continued fraction expansion of surd.

    Continued fraction expansion of a rational number is an expression
    obtained through an iterative process of representing a number as
    the sum of its integer part and the reciprocal of another number,
    then writing this other number as the sum of its integer part and
    another reciprocal, and so on.

    Here the three parameters are,
    * num: the rational part of the number's numerator
    * delta: the irrational part(discriminator) of the number's numerator
    * denominator: the denominator of the number

    ex: the golden ratio is (1 + sqrt(5))/2
    num = 1
    delta = 5
    den = 2

    The denominator of a rational number cannot be zero. So such
    input will result an error.

    >>> from sympy.ntheory.continued_fraction import continued_fraction
    >>> continued_fraction(1,0,0)
    Traceback (most recent call last):
    ...
    ValueError: The denominator is zero.

    If the discriminator is zero then the number will be a rational number.

    >>> from sympy.ntheory.continued_fraction import continued_fraction
    >>> continued_fraction(4,3,0)
    [1, 3]

    Golden ratio has the simplest continued fraction expansion,

    >>> from sympy.ntheory.continued_fraction import continued_fraction
    >>> continued_fraction(1,2,5)
    [[1], [1]]

    Note: if the length of recurrsive part of the continued part of
    expansion exceeds 100000 this module will truncate the convergents.

    See Also
    ========

    continued_fraction_rational_number
    sympy.ntheory.continued_fraction.continued_fraction_rational_number :
    function which calculate the continued fraction of a rational number.

    References
    ==========

    [1] A. J. van der Poorten, "NOTES ON CONTINUED FRACTIONS AND RECURRENCE
    SEQUENCES" in Number Theory and Cryptography. New York,
    USA: Cambridge university press, 2011,
    ch. 06, pp. 86-96.
    [2] http://en.wikipedia.org/wiki/Continued_fraction
    [3] http://www.numbertheory.org/ntw/N4.html#continued_fractions
    [4] http://www.numbertheory.org/pdfs/CFquadratic.pdf
    [5] http://maths.mq.edu.au/~alf/www-centre/alfpapers/a117.pdf
    """

    list = []

    # if the denominator is zero the expression cannot be a legal fraction
    if den == 0:
        raise ValueError("The denominator is zero.")

    # if the discriminator is negative the number is a complex number
    if delta < 0:
        raise ValueError("The number is not real, so it does not\n\
            have continued fraction expansion.")

    # if the discriminator is zero the number is a rational number
    if delta == 0:
        list = continued_fraction_rational_number(
            num, den)
        return list

    if (delta - (num*num)) % den != 0:
        delta = delta*den*den
        num = num*abs(den)
        den = den*abs(den)

    sqrtDelta = 0 | 1 << (delta.bit_length() + 1)/2
    sqrtDeltaNext = ((delta/sqrtDelta) + sqrtDelta) >> 1

    while sqrtDelta > sqrtDeltaNext:
        sqrtDelta = sqrtDeltaNext
        sqrtDeltaNext = ((delta/sqrtDelta) + sqrtDelta) >> 1

    if sqrtDelta*sqrtDelta == delta:
        continued_fraction_rational_number(num + sqrtDelta, den)
        return list

    biP = den

    if biP > 0:
        biK = sqrtDelta
    else:
        biK = sqrtDelta + 1

    biK = biK + num

    if biK > 0:
        if den > 0:
            biM = biK/den
        else:
            biM = ((den + 1) - biK)/(den*-1)
    else:
        if den > 0:
            biM = ((biK + 1) - den)/den
        else:
            biM = (biK*-1)/(den*-1)
    # appends the integer part of the continued fraction expansion
    # to the result list
    list.append(biM)
    biM = ((biM*den) - num)
    cont = -1
    K = -1
    P = -1
    L = -1
    M = -1

    while (cont < 0 or K != P or L != M):

        if (cont < 0 and biP > 0 and biP <= sqrtDelta + biM and
                biM > 0 and biM <= sqrtDelta):
            K = P = biP
            L = M = biM
            cont = 0

        # both numerator and denominator are positive
        if cont >= 0:
            P = (delta - (M*M))/P
            Z = (sqrtDelta + M)/P
            M = (Z*P) - M
            cont += 1
        else:
            biP = (delta - (biM*biM))/biP
            if biP > 0:
                Z = (sqrtDelta + biM)/biP
            else:
                Z = ((sqrtDelta + 1) + biM)/biP
            biM = (Z*biP) - biM

        # show convergent
        list.append(Z)

        # too many convergent
        if cont > 100000:
            raise NotImplementedError('more than 100,000 convergents')

    if cont >= 1:

        non_recurring_part = list[:len(list) - cont]

        recurring_part = list[len(list) - cont:]

        return [non_recurring_part, recurring_part]

    return list


def continued_fraction_rational_number(n, d):
    """This applies the continued fraction expansion to two numbers
    numerator/denominator

    >>> from sympy.ntheory.continued_fraction import\
    continued_fraction_rational_number
    >>> continued_fraction_rational_number(3, 8)
    [0, 2, 1, 2]
    """
    x = Rational(n, d)
    list = []
    while True:
        i = Integer(x)
        list.append(i)
        x -= i
        if not x:
            break
        x = 1/x
    return list
