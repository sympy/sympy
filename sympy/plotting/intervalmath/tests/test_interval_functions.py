from sympy.plotting.intervalmath import *


def test_functions():
    try:
        import numpy as np
    except ImportError:
        return
    else:
        sin_test()
        exp_test()
        log_test()
        atan_test()
        cos_test()
        tan_test()
        sqrt_test()
        imin_test()
        imax_test()
        cosh_test()
        sinh_test()
        tanh_test()
        asin_test()
        acos_test()
        floor_test()
        ceil_test()
        acos_test()
        asin_test()
        atanh_test()


def exp_test():
    import numpy as np
    a = exp(interval(-np.inf, 0))
    assert a.start == np.exp(-np.inf)
    assert a.end == np.exp(0)
    a = exp(interval(1, 2))
    assert a.start == np.exp(1)
    assert a.end == np.exp(2)


def log_test():
    import numpy as np
    a = log(interval(1, 2))
    assert a.start == 0
    assert a.end == np.log(2)
    a = log(interval(-1, 1))
    assert a.is_valid is None
    a = log(interval(-3, -1))
    assert a.is_valid is False


def atan_test():
    import numpy as np
    a = atan(interval(0, 1))
    assert a.start == np.arctan(0)
    assert a.end == np.arctan(1)
    a = atan(1)
    a.start == np.arctan(1)
    a.end == np.arctan(1)


def sin_test():
    import numpy as np
    a = sin(interval(0, np.pi / 4))
    assert a.start == np.sin(0)
    assert a.end == np.sin(np.pi / 4)

    a = sin(interval(-np.pi / 4, np.pi / 4))
    assert a.start == np.sin(-np.pi / 4)
    assert a.end == np.sin(np.pi / 4)

    a = sin(interval(np.pi / 4, 3 * np.pi / 4))
    assert a.start == np.sin(np.pi / 4)
    assert a.end == 1

    a = sin(interval(7 * np.pi / 6, 7 * np.pi / 4))
    assert a.start == -1
    assert a.end == np.sin(7 * np.pi / 6)

    a = sin(interval(0, 3 * np.pi))
    assert a.start == -1
    assert a.end == 1

    a = sin(interval(np.pi / 3, 7 * np.pi / 4))
    assert a.start == -1
    assert a.end == 1


def cos_test():
    import numpy as np
    a = cos(interval(0, np.pi / 4))
    assert a.start == np.cos(np.pi / 4)
    assert a.end == 1

    a = cos(interval(-np.pi / 4, np.pi / 4))
    assert a.start == np.cos(-np.pi / 4)
    assert a.end == 1

    a = cos(interval(np.pi / 4, 3 * np.pi / 4))
    assert a.start == np.cos(3 * np.pi / 4)
    assert a.end == np.cos(np.pi / 4)

    a = cos(interval(3 * np.pi / 4, 5 * np.pi / 4))
    assert a.start == -1
    assert a.end == np.cos(3 * np.pi / 4)

    a = cos(interval(0, 3 * np.pi))
    assert a.start == -1
    assert a.end == 1

    a = cos(interval(- np.pi / 3, 5 * np.pi / 4))
    assert a.start == -1
    assert a.end == 1


def tan_test():
    import numpy as np
    a = tan(interval(0, np.pi / 4))
    assert a.start == 0
    assert a.end == np.tan(np.pi / 4)

    a = tan(interval(np.pi / 4, 3 * np.pi / 4))
    #discontinuity
    assert a.is_valid is None


def sqrt_test():
    import numpy as np
    a = sqrt(interval(1, 4))
    assert a.start == 1
    assert a.end == 2

    a = sqrt(interval(0.01, 1))
    assert a.start == np.sqrt(0.01)
    assert a.end == 1

    a = sqrt(interval(-1, 1))
    assert a.is_valid is None

    a = sqrt(interval(-3, -1))
    assert a.is_valid is False


def imin_test():
    a = imin(interval(1, 3), interval(2, 5), interval(-1, 3))
    assert a.start == -1
    assert a.end == 3

    a = imin(-2, interval(1, 4))
    assert a.start == -2
    assert a.end == -2

    a = imin(5, interval(3, 4), interval(-2, 2, is_valid=False))
    assert a.start == 3
    assert a.end == 4


def imax_test():
    a = imax(interval(-2, 2), interval(2, 7), interval(-3, 9))
    assert a.start == 2
    assert a.end == 9

    a = imax(8, interval(1, 4))
    assert a.start == 8
    assert a.end == 8

    a = imax(interval(1, 2), interval(3, 4), interval(-2, 2, is_valid=False))
    assert a.start == 3
    assert a.end == 4


def sinh_test():
    import numpy as np
    a = sinh(interval(-1, 1))
    assert a.start == np.sinh(-1)
    assert a.end == np.sinh(1)


def cosh_test():
    import numpy as np
    a = cosh(interval(1, 2))
    assert a.start == np.cosh(1)
    assert a.end == np.cosh(2)
    a = cosh(interval(-2, -1))
    assert a.start == np.cosh(-1)
    assert a.end == np.cosh(-2)

    a = cosh(interval(-2, 1))
    assert a.start == 1
    assert a.end == np.cosh(-2)


def tanh_test():
    import numpy as np
    a = tanh(interval(-3, 3))
    assert a.start == np.tanh(-3)
    assert a.end == np.tanh(3)


def asin_test():
    import numpy as np
    a = asin(interval(-0.5, 0.5))
    assert a.start == np.arcsin(-0.5)
    assert a.end == np.arcsin(0.5)

    a = asin(interval(-1.5, 1.5))
    assert a.is_valid is None
    a = asin(interval(-2, -1.5))
    assert a.is_valid is False

    a = asin(interval(0, 2))
    assert a.is_valid is None

    a = asin(interval(2, 5))
    assert a.is_valid is False


def acos_test():
    import numpy as np
    a = acos(interval(-0.5, 0.5))
    assert a.start == np.arccos(0.5)
    assert a.end == np.arccos(-0.5)

    a = acos(interval(-1.5, 1.5))
    assert a.is_valid is None
    a = acos(interval(-2, -1.5))
    assert a.is_valid is False

    a = acos(interval(0, 2))
    assert a.is_valid is None

    a = acos(interval(2, 5))
    assert a.is_valid is False


def ceil_test():
    a = ceil(interval(0.2, 0.5))
    assert a.start == 1
    assert a.end == 1

    a = ceil(interval(0.5, 1.5))
    assert a.start == 1
    assert a.end == 2
    assert a.is_valid is None

    a = ceil(interval(-5, 5))
    assert a.is_valid is None


def floor_test():
    a = floor(interval(0.2, 0.5))
    assert a.start == 0
    assert a.end == 0

    a = floor(interval(0.5, 1.5))
    assert a.start == 0
    assert a.end == 1
    assert a.is_valid is None

    a = floor(interval(-5, 5))
    assert a.is_valid is None


def asinh_test():
    a = asinh(interval(1, 2))
    assert a.start == np.arcsinh(1)
    assert a.end == np.arcsinh(2)


def acosh_test():
    a = acosh(interval(3, 5))
    assert a.start == np.arccosh(3)
    assert a.end == np.arccosh(5)

    a = acosh(interval(0, 3))
    assert a.is_valid is None
    a = acosh(interval(-3, 0.5))
    assert a.is_valid is False


def atanh_test():
    a = atanh(interval(-0.5, 0.5))
    assert a.start == np.arctanh(-0.5)
    assert a.end == np.arctanh(0.5)

    a = atanh(interval(0, 3))
    assert a.is_valid is None

    a = atanh(interval(-3, -2))
    assert a.is_valid is False
