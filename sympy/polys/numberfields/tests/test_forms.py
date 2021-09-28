from sympy import QQ
from sympy.abc import x
from sympy.polys import Poly
from sympy.polys.numberfields.forms import StandardRep


def test_StandardRep_from_poly():
    T = Poly(x ** 2 + 5)
    f = Poly(QQ(3)*x/10 - QQ(7)/6)
    assert StandardRep.from_poly(T, f) == StandardRep(T, [-35, 9], 30)


def test_StandardRep_eq_int():
    T = Poly(x ** 2 + 5)
    assert StandardRep.zero(T) == 0
    assert StandardRep.one(T) == 1
    assert StandardRep(T, [2, 0], 1) == 2
    assert StandardRep(T, [10, 0], 5) == 2


def test_StandardRep_add():
    T = Poly(x ** 2 + 5)
    r = StandardRep(T, [1, 1], 6)
    s = StandardRep(T, [31, 1], 6)
    t = StandardRep(T, [5, 8], 10)
    u = StandardRep(T, [20, 29], 30)
    assert r + 5 == 5 + r == s
    assert r + t == t + r == u


def test_StandardRep_mul():
    T = Poly(x ** 2 + 5)
    r = StandardRep(T, [1, -2], 6)
    s = StandardRep(T, [-5, 10], 3)
    t = StandardRep(T, [95, 20], 18)
    assert -10 * r == r * (-10) == s
    assert r * s == s * r == t


def test_StandardRep_div():
    T = Poly(x ** 2 + 5)
    r = StandardRep(T, [4, -2], 6)
    s = StandardRep(T, [-2, 1], 18)
    assert r // -6 == s


def test_StandardRep_pow():
    T = Poly(x ** 2 + 5)
    r = StandardRep(T, [1, -2], 6)
    s = StandardRep(T, [-19, -4], 36)
    assert r ** 2 == s


def test_StandardRep_inverse():
    T = Poly(x ** 2 + 5)
    r = StandardRep(T, [1, -2], 6)
    assert r * (1 // r) == 1
    assert r * r**(-1) == 1


def test_StandardRep_mod_int():
    T = Poly(x ** 2 + 5)
    r = StandardRep(T, [1, 15], 2)
    s = StandardRep(T, [1, 1], 2)
    t = StandardRep(T, [1, 8], 2)
    assert r % 7 == s
    assert t % 7 == t
