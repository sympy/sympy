'''
This implementation is a heavily modified fixed point implementation of
BBP_formula for calculating the nth position of pi. The original hosted
at: http://en.literateprograms.org/Pi_with_the_BBP_formula_(Python)

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sub-license, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Modifications:
1.Once the nth digit is selected the number of digits of working
precision is calculated to ensure that the 14 Hexadecimal representation
of that region is accurate. This was found empirically to be
int((math.log10(n//1000))+18). This was found by searching for a value
of working precision for the n = 0 and n = 1 then n was increased until
the result was less precise, therefore increased again this was repeated
for increasing n and an effective fit was found between n and the
working precision value.

2. The while loop to evaluate whether the series has converged has be
replaced with a fixed for-loop, that option was selected because in a
very large number of cases the loop converged to a point where no
difference can be detected in less than 15 iterations. (done for more
accurate memory and time banking).

3. output hex string constrained to 14 characters (accuracy assured to be
n = 10**7)

4. pi_hex_digits(n) changed to have coefficient to the formula in an
array (perhaps just a matter of preference).

'''
from __future__ import print_function, division

import math
from sympy.core.compatibility import range


def _series(j, n, prec=14):

    # Left sum from the bbp algorithm
    s = 0
    D = _dn(n, prec)
    D4 = 4 * D
    k = 0
    d = 8 * k + j
    for k in range(n + 1):
        s += (pow(16, n - k, d) << D4) // d
        d += 8

    # Right sum iterates to infinity for full precision, but we
    # stop at the point where one iteration is beyond the precision
    # specified.

    t = 0
    k = n + 1
    e = 4*(D + n - k)
    d = 8 * k + j
    while True:
        dt = (1 << e) // d
        if not dt:
            break
        t += dt
        # k += 1
        e -= 4
        d += 8
    total = s + t

    return total


def pi_hex_digits(n, prec=14):
    """Returns a string containing ``prec`` (default 14) digits
    starting at the nth digit of pi in hex. Counting of digits
    starts at 0 and the decimal is not counted, so for n = 0 the
    returned value starts with 3; n = 1 corresponds to the first
    digit past the decimal point (which in hex is 2).

    Examples
    ========

    >>> from sympy.ntheory.bbp_pi import pi_hex_digits
    >>> pi_hex_digits(0)
    '3243f6a8885a30'
    >>> pi_hex_digits(13)
    '08d313198a2e03'

    References
    ==========

    .. [1] http://www.numberworld.org/digits/Pi/
    """

    # main of implementation arrays holding formulae coefficients
    n -= 1
    a = [4, 2, 1, 1]
    j = [1, 4, 5, 6]

    #formulae
    D = _dn(n, prec)
    x = + (a[0]*_series(j[0], n, prec)
         - a[1]*_series(j[1], n, prec)
         - a[2]*_series(j[2], n, prec)
         - a[3]*_series(j[3], n, prec)) & (16**D - 1)

    s = ("%0" + "%ix" % prec) % (x // 16**(D - prec))
    return s


def _dn(n, prec):
    # controller for n dependence on precision
    return int(math.log(n + prec)/math.log(16) + prec + 2)
