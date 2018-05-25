from __future__ import print_function, division

from sympy import sqrt, pi, E, exp
from sympy.core import S, Symbol, I
from sympy.core.compatibility import range
from sympy.discrete import convolution_fft, convolution_ntt
from sympy.utilities.pytest import raises
from sympy.abc import x, y

def test_convolution_fft():
    assert convolution_fft([1, 2, 3], [4, 5, 6]) == [4, 13, 28, 27, 18]
    assert convolution_fft([1], [5, 6, 7]) == [5, 6, 7]
    assert convolution_fft([1, 3], [5, 6, 7]) == [5, 21, 25, 21]

    assert convolution_fft([1 + 2*I], [2 + 3*I]) == [-4 + 7*I]
    assert convolution_fft([1 + 2*I, 3 + 4*I, 5 + S(3)/5*I], [S(2)/5 + S(4)/7*I]) == \
            [-S(26)/35 + 48*I/35, -S(38)/35 + 116*I/35, S(58)/35 + 542*I/175]

    assert convolution_fft([S(3)/4, S(5)/6], [S(7)/8, S(1)/3, S(2)/5]) == [S(21)/32,
                                                S(47)/48, S(26)/45, S(1)/3]
    assert convolution_fft([S(1)/9, S(2)/3, S(3)/5], [S(2)/5, S(3)/7, S(4)/9]) == [S(2)/45,
                                    S(11)/35, S(8152)/14175, S(523)/945, S(4)/15]
    assert convolution_fft([pi, E, sqrt(2)], [sqrt(3), 1/pi, 1/E]) == [sqrt(3)*pi,
                                                            1 + sqrt(3)*E,
                                                            E/pi + pi*exp(-1) + sqrt(6),
                                                            sqrt(2)/pi + 1,
                                                            sqrt(2)*exp(-1)]

    assert convolution_fft([2321, 33123], [5321, 6321, 71323]) == [12350041, 190918524,
                                                        374911166, 2362431729]
    assert convolution_fft([312313, 31278232], [32139631, 319631]) == [10037624576503,
                                                1005370659728895, 9997492572392]

    raises(TypeError, lambda: convolution_fft(x, y))
    raises(ValueError, lambda: convolution_fft([x, y], [y, x]))


def test_convolution_ntt():
    # prime moduli of the form (m*2**k + 1), sequence length
    # should be a divisor of 2**k
    p = 7*17*2**23 + 1
    q = 19*2**10 + 1
    r = 2*500000003 + 1 # only for sequences of length 1 or 2
    s = 2*3*5*7 # composite modulus

    assert convolution_ntt([2], [3], r) == [6]
    assert convolution_ntt([2, 3], [4], r) == [8, 12]

    assert convolution_ntt([32121, 42144, 4214, 4241], [32132, 3232, 87242], p) == [33867619,
                                    459741727, 79180879, 831885249, 381344700, 369993322]
    assert convolution_ntt([121913, 3171831, 31888131, 12], [17882, 21292, 29921, 312], q) == \
                                                [8158, 3065, 3682, 7090, 1239, 2232, 3744]

    assert convolution_ntt([12, 19, 21, 98, 67], [2, 6, 7, 8, 9], p) == \
                    convolution_ntt([12, 19, 21, 98, 67], [2, 6, 7, 8, 9], q)
    assert convolution_ntt([12, 19, 21, 98, 67], [21, 76, 17, 78, 69], p) == \
                    convolution_ntt([12, 19, 21, 98, 67], [21, 76, 17, 78, 69], q)

    raises(ValueError, lambda: convolution_ntt([2, 3], [4, 5], r))
    raises(ValueError, lambda: convolution_ntt([x, y], [y, x], q))
    raises(TypeError, lambda: convolution_ntt(x, y, p))
