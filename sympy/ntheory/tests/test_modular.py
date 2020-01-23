from sympy.ntheory.modular import crt, crt1, crt2, solve_congruence, \
    crt_cartesian
from sympy.testing.pytest import raises
from sympy.utilities.iterables import cartes


def test_crt():
    def mcrt(m, v, r, symmetric=False):
        assert crt(m, v, symmetric)[0] == r
        mm, e, s = crt1(m)
        assert crt2(m, v, mm, e, s, symmetric) == (r, mm)

    mcrt([2, 3, 5], [0, 0, 0], 0)
    mcrt([2, 3, 5], [1, 1, 1], 1)

    mcrt([2, 3, 5], [-1, -1, -1], -1, True)
    mcrt([2, 3, 5], [-1, -1, -1], 2*3*5 - 1, False)

    assert crt([656, 350], [811, 133], symmetric=True) == (
        -56917, 114800)

    raises(ValueError, lambda: crt(
        [20, 2, 5], [1, 1, 2]))


def test_modular():
    assert solve_congruence(*list(zip([3, 4, 2], [12, 35, 17]))
        ) == (1719, 7140)
    assert solve_congruence(*list(zip([3, 4, 2], [12, 6, 17]))) is None
    assert solve_congruence(*list(zip([3, 4, 2], [13, 7, 17]))
        ) == (172, 1547)
    assert solve_congruence(*list(zip([-10, -3, -15], [13, 7, 17]))
        ) == (172, 1547)
    assert solve_congruence(*list(
        zip([-10, -3, 1, -15], [13, 7, 7, 17]))) is None
    assert solve_congruence(
        *list(zip([-10, -5, 2, -15], [13, 7, 7, 17]))) == (835, 1547)
    assert solve_congruence(
        *list(zip([-10, -5, 2, -15], [13, 7, 14, 17]))) == (2382, 3094)
    assert solve_congruence(
        *list(zip([-10, 2, 2, -15], [13, 7, 14, 17]))) == (2382, 3094)
    assert solve_congruence(*list(zip((1, 1, 2), (3, 2, 4)))) is None
    raises(
        ValueError, lambda: solve_congruence(*list(zip(
        [3, 4, 2], [12.1, 35, 17]))))


def test_crt_cartesian():
    def check(r, p, ans):
        a = crt_cartesian(r, p)
        assert a == ans
        for i in zip([tuple([n%i for i in p]) for n in a], cartes(*r)):
            print(i)
            assert i[0] == i[1]

    check([[3, 5], [3, 7]], [7, 11], [3, 73, 47, 40])
    check([[3], [3, 7]], [7, 11], [3, 73])
    check([[1, 5], [4, 7], [6, 8]], [6, 11, 13],
        [565, 697, 799, 73, 851, 125, 227, 359])
    check([[5, 3], [7, 9], [3, 7]], [7, 11, 13],
        [887, 579, 614, 306, 458, 150, 185, 878])
    check([[11, 51], [54, 72], [16, 38]], [67, 79, 43],
        [196053, 111365, 129790, 45102, 189259, 104571, 122996, 38308])
    raises(ValueError, lambda : crt_cartesian(
        [[4, 2], [5, 7], [1, 3]], [6, 11, 14]))
    raises(ValueError, lambda: crt_cartesian(
        [[2, 3], [3, 4]], [4, 5, 6]))
    raises(ValueError, lambda: crt_cartesian([[4, 7], [3, 5]], []))
    raises(ValueError, lambda: crt_cartesian([[4, 7], []], [9, 11]))
