from sympy.external.ntheory import bit_scan1, bit_scan0


def test_bit_scan1():
    assert bit_scan1(0) is None
    assert bit_scan1(1) == 0
    assert bit_scan1(-1) == 0
    assert bit_scan1(2) == 1
    assert bit_scan1(7) == 0
    assert bit_scan1(-7) == 0
    for i in range(100):
        assert bit_scan1(1 << i) == i
        assert bit_scan1((1 << i) * 31337) == i
    for i in range(500):
        n = (1 << 500) + (1 << i)
        assert bit_scan1(n) == i
    assert bit_scan1(1 << 1000001) == 1000001
    assert bit_scan1((1 << 273956)*7**37) == 273956
    # issue 12709
    for i in range(1, 10):
        big = 1 << i
        assert bit_scan1(-big) == bit_scan1(big)


def test_bit_scan0():
    assert bit_scan0(-1) is None
    assert bit_scan0(0) == 0
    assert bit_scan0(1) == 1
    assert bit_scan0(-2) == 0
