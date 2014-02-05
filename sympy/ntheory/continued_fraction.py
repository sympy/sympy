from fractions import gcd


def continued_fraction(numerator, denominator, delta):
    """
        Return continued fraction expansion of surd.

        Continued fraction expansion of a rational number is an expression
        obtained through an iterative process of representing a number as
        the sum of its integer part and the reciprocal of another number,
        then writing this other number as the sum of its integer part and
        another reciprocal, and so on.

        Here the three parameters are,
        * numerator: the rational part of the number's numerator
        * delta: the irrational part(discriminator) of the number's numerator
        * denominator: the denominator of the number

        ex: the golden ratio is (1+sqrt(5))/2
        numerator = 1
        delta = 5
        denominator = 2

        The denominator of a rational number cannot be zero. So such
        input will result an error.

        >>> continued_fraction(1,0,0)
        'Error: The denominator is zero.'

        If the discriminator is negative the number is a complex number.
        Complex numbers does not have continued fraction expansions.

        >>> continued_fraction(1,2,-1)
        'Error: The number is not real, so it does not have continued fraction
        expansion.'

        If the discriminator is zero then the number will be a rational number.

        >>> continued_fraction(4,3,0)
        [1, 3]

        Golden ratio has the simplest continued fraction expansion,

        >>> continued_fraction(1,2,5)
        [[1], [1]]

        Note: if the length of recurrsive part of the continued part of
        expantion exceeds 100000 this module will truncate the convergents.

        See Also
        ========

        continued_fraction_rational_number
        sympy.ntheory.continued_fraction.continued_fraction_rational_number :
        function which calculate the continued fraction of a rational number.

        References
        ==========

        - A. J. van der Poorten, "NOTES ON CONTINUED FRACTIONS AND RECURRENCE
        SEQUENCES" in Number Theory and Cryptography. New York,
        USA: Cambridge university press, 2011,
        ch. 06, pp. 86-96.
        - http://en.wikipedia.org/wiki/Continued_fraction
        - http://www.numbertheory.org/ntw/N4.html#continued_fractions
        - http://www.numbertheory.org/pdfs/CFquadratic.pdf
        - http://maths.mq.edu.au/~alf/www-centre/alfpapers/a117.pdf
    """

    _continued_fraction_expansion = []

    #if the denominator is zero the expression cannot be a legal fraction
    if denominator == 0:
        raise ValueError("The denominator is zero.")

    #if the discriminator is negative the number is a complex number
    if delta < 0:
        raise ValueError("The number is not real, so it does not\n\
            have continued fraction expansion.")

    #if the discriminator is zero the number is a rational number
    if delta == 0:
        _continued_fraction_expansion = continued_fraction_rational_number(
            numerator, denominator)
        return _continued_fraction_expansion

    if (delta-(numerator*numerator)) % denominator != 0:
        delta = delta*denominator*denominator
        numerator = numerator*abs(denominator)
        denominator = denominator*abs(denominator)

    sqrtDelta = 0 | 1 << (delta.bit_length()+1)/2
    sqrtDeltaNext = ((delta/sqrtDelta)+sqrtDelta) >> 1

    while sqrtDelta > sqrtDeltaNext:
        sqrtDelta = sqrtDeltaNext
        sqrtDeltaNext = ((delta/sqrtDelta)+sqrtDelta) >> 1

    if sqrtDelta*sqrtDelta == delta:
        continued_fraction_rational_number(numerator+sqrtDelta, denominator)
        return _continued_fraction_expansion

    biP = denominator

    if biP > 0:
        biK = sqrtDelta
    else:
        biK = sqrtDelta+1

    biK = biK+numerator

    if biK > 0:
        if denominator > 0:
            biM = biK/denominator
        else:
            biM = ((denominator+1)-biK)/(denominator*-1)
    else:
        if denominator > 0:
            biM = ((biK+1)-denominator)/denominator
        else:
            biM = (biK*-1)/(denominator*-1)
    #appends the integer part of the continued fraction expansion
    #to the result list
    _continued_fraction_expansion.append(biM)
    biM = ((biM*denominator)-numerator)
    cont = -1
    K = -1
    P = -1
    L = -1
    M = -1

    while (cont < 0 or K != P or L != M):

        if (cont < 0 and biP > 0 and biP <= sqrtDelta+biM and
                biM > 0 and biM <= sqrtDelta):
            K = P = biP
            L = M = biM
            cont = 0

        #both numerator and denominator are positive
        if cont >= 0:
            P = (delta-(M*M))/P
            Z = (sqrtDelta+M)/P
            M = (Z*P)-M
            cont += 1
        else:
            biP = (delta-(biM*biM))/biP
            if biP > 0:
                Z = (sqrtDelta+biM)/biP
            else:
                Z = ((sqrtDelta+1)+biM)/biP
            biM = (Z*biP)-biM

        #show convergent
        _continued_fraction_expansion.append(Z)

        #too many convergent
        if cont > 100000:
            return "truncate after 100000 convergents"

    if cont >= 1:

        non_recurring_part = _continued_fraction_expansion[
            :len(_continued_fraction_expansion)-cont]

        recurring_part = _continued_fraction_expansion[
            len(_continued_fraction_expansion)-cont:]

        return [non_recurring_part, recurring_part]

    return _continued_fraction_expansion


def continued_fraction_rational_number(numerator, denominator):
    """This applies the continued fraction expansion to two numbers
    numerator/denominator

    >>> continued_fraction(3, 8)
    [0, 2, 1, 2]
    """
    numerator = int(numerator)
    denominator = int(denominator)
    temp = numerator//denominator
    if temp*denominator == numerator:
        return [temp, ]

    _continued_fraction_expansion = continued_fraction_rational_number(
        denominator, numerator - temp*denominator)
    _continued_fraction_expansion.insert(0, temp)
    return _continued_fraction_expansion
