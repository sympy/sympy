from __future__ import print_function, division

from sympy import sqrt, pi, E, exp
from sympy.core import S, Symbol, I
from sympy.core.compatibility import range
from sympy.discrete import convolution, convolution_fft, convolution_ntt
from sympy.utilities.pytest import raises
from sympy.abc import x, y

def test_convolution():
    a = [1, S(5)/3, sqrt(3), S(7)/5]
    b = [9, 5, 5, 4, 3, 2]
    c = [3, 5, 3, 7, 8]
    d = [1422, 6572, 3213, 5552]

    assert convolution(a, b, fft=True) == convolution_fft(a, b)
    assert convolution(a, b, dps=9, fft=True) == convolution_fft(a, b, dps=9)
    assert convolution(a, d, fft=True, dps=7) == convolution_fft(d, a, dps=7)
    assert convolution(a, d[1:], dps=3) == convolution_fft(d[1:], a, dps=3)

    # prime moduli of the form (m*2**k + 1), sequence length
    # should be a divisor of 2**k
    p = 7*17*2**23 + 1
    q = 19*2**10 + 1
    assert convolution(d, b, ntt=True, prime=q) == convolution_ntt(b, d, prime=q)
    assert convolution(c, b, prime=p) == convolution_ntt(b, c, prime=p)
    assert convolution(d, c, prime=p, ntt=True) == convolution_ntt(c, d, prime=p)

    raises(TypeError, lambda: convolution(b, d, ntt=True))
    raises(TypeError, lambda: convolution(b, d, ntt=True, cycle=None))
    raises(TypeError, lambda: convolution(b, d, dps=5, prime=q))
    raises(TypeError, lambda: convolution(b, d, dps=6, ntt=True, prime=q))
    raises(TypeError, lambda: convolution(b, d, fft=True, dps=7, ntt=True, prime=q))

    # ntt is a specialized variant of fft, TypeError should not be raised
    assert convolution(b, d, fft=True, ntt=True, prime=q) == \
            convolution_ntt(b, d, prime=q)


def test_cyclic_convolution():
    # fft
    a = [1, S(5)/3, sqrt(3), S(7)/5]
    b = [9, 5, 5, 4, 3, 2]

    assert convolution([1, 2, 3], [4, 5, 6], cycle=None) == \
            convolution([1, 2, 3], [4, 5, 6], cycle=5) == \
                convolution([1, 2, 3], [4, 5, 6])

    assert convolution([1, 2, 3], [4, 5, 6], cycle=3) == [31, 31, 28]

    assert convolution(a, b, fft=True, cycle=4) == \
            convolution(a, b, cycle=4)

    assert convolution(a, b, fft=True, dps=3, cycle=4) == \
            convolution(a, b, dps=3, cycle=4)

    a = [S(1)/3, S(7)/3, S(5)/9, S(2)/7, S(5)/8]
    b = [S(3)/5, S(4)/7, S(7)/8, S(8)/9]

    assert convolution(a, b, cycle=None) == \
            convolution(a, b, cycle=len(a) + len(b) - 1)

    assert convolution(a, b, cycle=4) == [S(87277)/26460, S(30521)/11340,
                            S(11125)/4032, S(3653)/1080]

    assert convolution(a, b, cycle=6) == [S(20177)/20160, S(676)/315, S(47)/24,
                            S(3053)/1080, S(16397)/5292, S(2497)/2268]

    assert convolution(a, b, cycle=9) == \
                convolution(a, b, cycle=None) + [S.Zero]

    # ntt
    a = [2313, 5323532, S(3232), 42142, 42242421]
    b = [S(33456), 56757, 45754, 432423]

    assert convolution(a, b, prime=19*2**10 + 1) == \
            convolution(a, b, prime=19*2**10 + 1, cycle=None) == \
                convolution(a, b, prime=19*2**10 + 1, cycle=8)

    assert convolution(a, b, prime=19*2**10 + 1, cycle=5) == [96, 17146, 2664,
                                                                    15534, 3517]

    assert convolution(a, b, prime=19*2**10 + 1, cycle=7) == [4643, 3458, 1260,
                                                        15534, 3517, 16314, 13688]

    assert convolution(a, b, prime=19*2**10 + 1, cycle=9) == \
            convolution(a, b, prime=19*2**10 + 1) + [0]


def test_convolution_fft():
    assert all(convolution_fft([], x, dps=y) == [] for x in ([], [1]) for y in (None, 3))
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

    assert all(convolution_ntt([], x, prime=y) == [] for x in ([], [1]) for y in (p, q, r))
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
