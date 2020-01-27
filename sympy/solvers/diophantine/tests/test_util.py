from sympy import (symbols, factorint, oo)
from sympy.solvers.diophantine.diophantine import (diop_general_sum_of_squares)
from sympy.solvers.diophantine.util import (length, square_factor, partition, prime_as_sum_of_two_squares,
    sum_of_three_squares, sum_of_four_squares, power_representation, sum_of_squares, can_do_sum_of_squares,
    gaussian_reduce, nint_or_floor, odd, even, remove_gcd, descent, ldescent)
from sympy.testing.pytest import raises


def test_length():
    assert length(2, 1, 0) == 1
    assert length(-2, 4, 5) == 3
    assert length(-5, 4, 17) == 4
    assert length(0, 4, 13) == 6
    assert length(7, 13, 11) == 23
    assert length(1, 6, 4) == 2


def test_square_factor():
    assert square_factor(1) == square_factor(-1) == 1
    assert square_factor(0) == 1
    assert square_factor(5) == square_factor(-5) == 1
    assert square_factor(4) == square_factor(-4) == 2
    assert square_factor(12) == square_factor(-12) == 2
    assert square_factor(6) == 1
    assert square_factor(18) == 3
    assert square_factor(52) == 2
    assert square_factor(49) == 7
    assert square_factor(392) == 14
    assert square_factor(factorint(-12)) == 2


def test_diop_partition():
    for n in [8, 10]:
        for k in range(1, 8):
            for p in partition(n, k):
                assert len(p) == k
    assert [p for p in partition(3, 5)] == []
    assert [list(p) for p in partition(3, 5, 1)] == [
        [0, 0, 0, 0, 3], [0, 0, 0, 1, 2], [0, 0, 1, 1, 1]]
    assert list(partition(0)) == [()]
    assert list(partition(1, 0)) == [()]
    assert [list(i) for i in partition(3)] == [[1, 1, 1], [1, 2], [3]]


def test_prime_as_sum_of_two_squares():
    for i in [5, 13, 17, 29, 37, 41, 2341, 3557, 34841, 64601]:
        a, b = prime_as_sum_of_two_squares(i)
        assert a ** 2 + b ** 2 == i
    assert prime_as_sum_of_two_squares(7) is None
    ans = prime_as_sum_of_two_squares(800029)
    assert ans == (450, 773) and type(ans[0]) is int


def test_sum_of_three_squares():
    for i in [0, 1, 2, 34, 123, 34304595905, 34304595905394941, 343045959052344,
              800, 801, 802, 803, 804, 805, 806]:
        a, b, c = sum_of_three_squares(i)
        assert a ** 2 + b ** 2 + c ** 2 == i

    assert sum_of_three_squares(7) is None
    assert sum_of_three_squares((4 ** 5) * 15) is None
    assert sum_of_three_squares(25) == (5, 0, 0)
    assert sum_of_three_squares(4) == (0, 0, 2)


def test_sum_of_four_squares():
    from random import randint

    # this should never fail
    n = randint(1, 100000000000000)
    assert sum(i ** 2 for i in sum_of_four_squares(n)) == n

    assert sum_of_four_squares(0) == (0, 0, 0, 0)
    assert sum_of_four_squares(14) == (0, 1, 2, 3)
    assert sum_of_four_squares(15) == (1, 1, 2, 3)
    assert sum_of_four_squares(18) == (1, 2, 2, 3)
    assert sum_of_four_squares(19) == (0, 1, 3, 3)
    assert sum_of_four_squares(48) == (0, 4, 4, 4)


def test_power_representation():
    tests = [(1729, 3, 2), (234, 2, 4), (2, 1, 2), (3, 1, 3), (5, 2, 2), (12352, 2, 4),
             (32760, 2, 3)]

    for test in tests:
        n, p, k = test
        f = power_representation(n, p, k)

        while True:
            try:
                l = next(f)
                assert len(l) == k

                chk_sum = 0
                for l_i in l:
                    chk_sum = chk_sum + l_i ** p
                assert chk_sum == n

            except StopIteration:
                break

    assert list(power_representation(20, 2, 4, True)) == \
           [(1, 1, 3, 3), (0, 0, 2, 4)]
    raises(ValueError, lambda: list(power_representation(1.2, 2, 2)))
    raises(ValueError, lambda: list(power_representation(2, 0, 2)))
    raises(ValueError, lambda: list(power_representation(2, 2, 0)))
    assert list(power_representation(-1, 2, 2)) == []
    assert list(power_representation(1, 1, 1)) == [(1,)]
    assert list(power_representation(3, 2, 1)) == []
    assert list(power_representation(4, 2, 1)) == [(2,)]
    assert list(power_representation(3 ** 4, 4, 6, zeros=True)) == \
           [(1, 2, 2, 2, 2, 2), (0, 0, 0, 0, 0, 3)]
    assert list(power_representation(3 ** 4, 4, 5, zeros=False)) == []
    assert list(power_representation(-2, 3, 2)) == [(-1, -1)]
    assert list(power_representation(-2, 4, 2)) == []
    assert list(power_representation(0, 3, 2, True)) == [(0, 0)]
    assert list(power_representation(0, 3, 2, False)) == []
    # when we are dealing with squares, do feasibility checks
    assert len(list(power_representation(4 ** 10 * (8 * 10 + 7), 2, 3))) == 0
    # there will be a recursion error if these aren't recognized
    big = 2 ** 30
    for i in [13, 10, 7, 5, 4, 2, 1]:
        assert list(power_representation(big, 2, big - i)) == []


def test_sum_of_squares_powers():
    tru = set([
        (0, 0, 1, 1, 11), (0, 0, 5, 7, 7), (0, 1, 3, 7, 8), (0, 1, 4, 5, 9),
        (0, 3, 4, 7, 7), (0, 3, 5, 5, 8), (1, 1, 2, 6, 9), (1, 1, 6, 6, 7),
        (1, 2, 3, 3, 10), (1, 3, 4, 4, 9), (1, 5, 5, 6, 6), (2, 2, 3, 5, 9),
        (2, 3, 5, 6, 7), (3, 3, 4, 5, 8)])

    u, v, x, y, z = symbols('u, v, x, y, z', integer=True)
    eq = u ** 2 + v ** 2 + x ** 2 + y ** 2 + z ** 2 - 123
    ans = diop_general_sum_of_squares(eq, oo)  # allow oo to be used
    assert len(ans) == 14
    assert ans == tru

    raises(ValueError, lambda: list(sum_of_squares(10, -1)))
    assert list(sum_of_squares(-10, 2)) == []
    assert list(sum_of_squares(2, 3)) == []
    assert list(sum_of_squares(0, 3, True)) == [(0, 0, 0)]
    assert list(sum_of_squares(0, 3)) == []
    assert list(sum_of_squares(4, 1)) == [(2,)]
    assert list(sum_of_squares(5, 1)) == []
    assert list(sum_of_squares(50, 2)) == [(5, 5), (1, 7)]
    assert list(sum_of_squares(11, 5, True)) == [
        (1, 1, 1, 2, 2), (0, 0, 1, 1, 3)]
    assert list(sum_of_squares(8, 8)) == [(1, 1, 1, 1, 1, 1, 1, 1)]

    assert [len(list(sum_of_squares(i, 5, True))) for i in range(30)] == [
        1, 1, 1, 1, 2,
        2, 1, 1, 2, 2,
        2, 2, 2, 3, 2,
        1, 3, 3, 3, 3,
        4, 3, 3, 2, 2,
        4, 4, 4, 4, 5]
    assert [len(list(sum_of_squares(i, 5))) for i in range(30)] == [
        0, 0, 0, 0, 0,
        1, 0, 0, 1, 0,
        0, 1, 0, 1, 1,
        0, 1, 1, 0, 1,
        2, 1, 1, 1, 1,
        1, 1, 1, 1, 3]
    for i in range(30):
        s1 = set(sum_of_squares(i, 5, True))
        assert not s1 or all(sum(j ** 2 for j in t) == i for t in s1)
        s2 = set(sum_of_squares(i, 5))
        assert all(sum(j ** 2 for j in t) == i for t in s2)

    raises(ValueError, lambda: list(power_representation(2, -1, 1)))
    raises(ValueError, lambda: list(power_representation(2, 1, -1)))
    assert list(power_representation(-2, 3, 2)) == [(-1, -1)]
    assert list(power_representation(-2, 4, 2)) == []
    assert list(power_representation(2, 1, 1)) == [(2,)]
    assert list(power_representation(2, 1, 3, True)) == [(0, 0, 2), (0, 1, 1)]
    assert list(power_representation(5, 1, 2, True)) == [(0, 5), (1, 4), (2, 3)]
    assert list(power_representation(6, 2, 2)) == []
    assert list(power_representation(3 ** 5, 3, 1)) == []
    assert list(power_representation(3 ** 6, 3, 1)) == [(9,)] and (9 ** 3 == 3 ** 6)
    assert list(power_representation(2 ** 1000, 5, 2)) == []


def test__can_do_sum_of_squares():
    assert can_do_sum_of_squares(3, -1) is False
    assert can_do_sum_of_squares(-3, 1) is False
    assert can_do_sum_of_squares(0, 1)
    assert can_do_sum_of_squares(4, 1)
    assert can_do_sum_of_squares(1, 2)
    assert can_do_sum_of_squares(2, 2)
    assert can_do_sum_of_squares(3, 2) is False


def test_gaussian_reduce():
    assert gaussian_reduce(4, 1, 3) == (1, 1)


def test_nint_or_floor():
    assert nint_or_floor(16, 10) == 2


def test_odd_even():
    assert odd(1) == (not even(1)) == True
    assert odd(0) == (not even(0)) == False


def test_remove_gcd():
    assert remove_gcd(2, 4, 6) == (1, 2, 3)
    raises(TypeError, lambda: remove_gcd((2, 4, 6)))


def test_descent():
    u = ([(13, 23), (3, -11), (41, -113), (91, -3), (1, 1), (1, -1), (17, 13), (123689, 1), (19, -570)])
    for a, b in u:
        w, x, y = descent(a, b)
        assert a*x**2 + b*y**2 == w**2
    # the docstring warns against bad input, so these are expected results
    # - can't both be negative
    raises(TypeError, lambda: descent(-1, -3))
    # A can't be zero unless B != 1
    raises(ZeroDivisionError, lambda: descent(0, 3))
    # supposed to be square-free
    raises(TypeError, lambda: descent(4, 3))


def test_ldescent():
    # Equations which have solutions
    u = ([(13, 23), (3, -11), (41, -113), (4, -7), (-7, 4), (91, -3), (1, 1), (1, -1),
        (4, 32), (17, 13), (123689, 1), (19, -570)])
    for a, b in u:
        w, x, y = ldescent(a, b)
        assert a*x**2 + b*y**2 == w**2
    assert ldescent(-1, -1) is None
