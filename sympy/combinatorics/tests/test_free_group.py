from sympy.combinatorics.free_group import free_group, FreeGroup
from sympy.core import Symbol
from sympy.utilities.pytest import raises
from sympy import oo

F, x, y, z = free_group("x, y, z")


def test_FreeGroup__init__():
    x, y, z = map(Symbol, "xyz")

    assert len(FreeGroup("x, y, z").generators) == 3
    assert len(FreeGroup(x).generators) == 1
    assert len(FreeGroup(("x", "y", "z"))) == 3
    assert len(FreeGroup((x, y, z)).generators) == 3


def test_free_group():
    G, a, b, c = free_group("a, b, c")
    assert F.generators == (x, y, z)
    assert x*z**2 in F
    assert x in F
    assert y*z**-1 in F
    assert (y*z)**0 in F
    assert a not in F
    assert a**0 not in F
    assert len(F) == 3
    assert str(F) == '<free group on the generators (x, y, z)>'
    assert not F == G
    assert F.order() == oo
    assert F.is_abelian == False
    assert F.center() == set([F.identity])

    (e,) = free_group("")
    assert e.order() == 1
    assert e.generators == ()
    assert e.elements == set([e.identity])
    assert e.is_abelian == True


def test_FreeGroup__hash__():
    assert hash(F)


def test_FreeGroup__eq__():
    assert free_group("x, y, z")[0] == free_group("x, y, z")[0]
    assert free_group("x, y, z")[0] is free_group("x, y, z")[0]

    assert free_group("x, y, z")[0] != free_group("a, x, y")[0]
    assert free_group("x, y, z")[0] is not free_group("a, x, y")[0]

    assert free_group("x, y")[0] != free_group("x, y, z")[0]
    assert free_group("x, y")[0] is not free_group("x, y, z")[0]

    assert free_group("x, y, z")[0] != free_group("x, y")[0]
    assert free_group("x, y, z")[0] is not free_group("x, y")[0]


def test_FreeGroup__getitem__():
    assert F[0:] == FreeGroup("x, y, z")
    assert F[1:] == FreeGroup("y, z")
    assert F[2:] == FreeGroup("z")
    # modify
    #assert F[3:] == FreeGroup("")


def test_FreeGroupElm__hash__():
    assert hash(x*y*z)


def test_FreeGroupElm_copy():
    f = x*y*z**3
    g = f.copy()
    h = x*y*z**7

    assert f == g
    assert f != h
