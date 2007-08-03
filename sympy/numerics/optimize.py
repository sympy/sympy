from sympy import Polynomial, Symbol
from float_ import Float, ComplexFloat
from evalf_ import polyfunc


def polyroots(poly, maxsteps=20):
    """
    Numerically locate all (complex) roots of a given polynomial 'poly'
    using the Durand-Kerner method (http://en.wikipedia.org/wiki/
    Durand-Kerner_method)

    This function returns a tuple (roots, err) where roots is a list of
    ComplexFloats sorted by absolute value, and err is an estimate of
    the maximum error. The 'poly' argument should be a SymPy expression
    representing a univariate polynomial.

    Example:
        >>> Float.setdps(15)
        >>> x = Symbol('x')
        >>> r, e = polyroots((x-3)*(x-2))
        >>> r[0]
        ComplexFloat(real='1.9999999999999993', imag='5.1698788284564230E-26')
        >>> r[1]
        ComplexFloat(real='3', imag='-4.6222318665293660E-33')
        >>> e
        Float('7.0435024728931298E-18')

    Polyroots attempts to achieve a maximum error less than the epsilon
    of the current working precision, but may fail to do so if the
    polynomial is badly conditioned. Usually the error can be reduced
    by increasing the 'maxsteps' parameter, although this will of
    course also reduce the speed proprtionally. It may also be
    necessary to increase the working precision.

    In particular, polyroots is easily fooled by multiple roots and
    roots that are very close together. In general, if n-multiple roots
    are present, polyroots will typically only locate their values to
    within 1/n of the working precision. For example, with the standard
    precision of ~15 decimal digits, the expected accuracy for a double
    root is 7-8 decimals:

        >>> Float.setdps(15)
        >>> r, e = polyroots((x-2)**2)
        >>> r[0]
        ComplexFloat(real='1.9999999832018203', imag='-2.5469994476319085E-8')
        >>> r[1]
        ComplexFloat(real='2.0000000110216827', imag='1.4784776969781909E-8')
        >>> e
        Float('5.4076851354766810E-8')

    The error estimate is also easily fooled when multiple roots are
    present. In the previous example, increasing 'maxsteps' will only
    have this effect:

        >>> r, e = polyroots((x-2)**2, maxsteps=30)
        >>> r[0].real
        Float('1.9999999832018203')
        >>> e
        Float('0')

    An effective cure is to multiply the working precision n times (in
    the case of a double root, doubling it):

        >>> Float.setdps(30)
        >>> r, e = polyroots((x-2)**2, 40)
        >>> r[0].real
        Float('1.9999999999999997562825649058053')
        >>> e
        Float('0')

    The error estimate is again wrong, but the result is good to 15
    decimals.
    """
    # Must be monic
    poly = Polynomial(poly)
    poly = Polynomial(poly / poly.coeffs[0][0])
    deg = poly.coeffs[0][1]
    f = polyfunc(poly)
    roots = [ComplexFloat(0.4+0.9j)**n for n in range(deg)]
    error = [Float(1) for n in range(deg)]
    for step in range(maxsteps):
        if max(error).ae(0):
            break
        for i in range(deg):
            if not error[i].ae(0):
                p = roots[i]
                x = f(p)
                for j in range(deg):
                    if i != j:
                        try:
                            x /= (p-roots[j])
                        except ZeroDivisionError:
                            continue
                roots[i] = p - x
                error[i] = abs(x)
    roots.sort(key=abs)
    return roots, max(error)
