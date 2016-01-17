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
replaced with a fixed for loop, that option was selected because in a
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


def _series(j, n):

    # Left sum from the bbp algorithm
    s = 0
    D = _dn(n)
    for k in range(n + 1):
        r = 8*k + j
        s += (pow(16, n - k, r) << 4 * D) // r

    # Right sum. should iterate to infinty, but now just iterates to the point where
    # one iterations change is beyond the resolution of the data type used

    t = 0
    for k in range(n + 1, n + 15):
        xp = int(16**(n - k) * 16**D)
        t += xp // (8 * k + j)
    total = s + t

    return total


def pi_hex_digits(n):
    """Returns a string containing 14 digits after the nth value of pi in hex
       The decimal has been taken out of the number, so
       n = 0[0] = 3 # First digit of pi in hex, 3

    Examples
    ========

    >>> from sympy.ntheory.bbp_pi import pi_hex_digits
    >>> pi_hex_digits(0)
    '3243f6a8885a30'
    >>> pi_hex_digits(10)
    '5a308d313198a2'
    """

    # main of implementation arrays holding formulae coefficients
    n -= 1
    a = [4, 2, 1, 1]
    j = [1, 4, 5, 6]

    #formulae
    x = + (a[0]*_series(j[0], n)
         - a[1]*_series(j[1], n)
         - a[2]*_series(j[2], n)
         - a[3]*_series(j[3], n)) & (16**(_dn(n)) - 1)

    s = ("%016x" % x)
    #s is constrained between 0 and 14
    return s[:14]


def _dn(n):
    # controller for n dependence on precision
    if (n < 1000):
        return 16
    return int(math.log10(n//1000) + 18)
