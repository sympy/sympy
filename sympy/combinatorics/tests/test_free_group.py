from sympy.combinatorics.free_group import FreeGroup
from sympy.utilities.pytest import raises
from sympy import S

f = FreeGroup(4)

def test_method():
    assert f.generators == [f[0], f[1], f[2], f[3]]
    assert len(f) == 4
    raises(IndexError, lambda: f[4])
    assert str(f) == '<free group on the generators [f0, f1, f2, f3]>'
    assert f.order() == S.Infinity


def test_contains():
    g = FreeGroup(4)
    assert f[0] in f
    assert not f[0] in g


def test_freegroupelm():
    assert (f[0]**3).letter_form == list([1, 1, 1])
    assert (f[0]**-2*f[1]**3*f[3]).letter_form == list([-1, -1, 2, 2, 2, 4])
    assert ((f[1])**0).order() == 1
