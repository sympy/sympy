from sympy import symbols
from sympy.external import import_module
from sympy.utilities.matchpy_connector import WildDot, WildPlus, WildStar

matchpy = import_module("matchpy")

x, y, z = symbols("x y z")


def _get_first_match(expr, pattern):
    from matchpy import ManyToOneMatcher, Pattern

    matcher = ManyToOneMatcher()
    matcher.add(Pattern(pattern))
    return next(iter(matcher.match(expr)))


def test_matchpy_connector():
    if matchpy is None:
        return

    from multiset import Multiset
    from matchpy import Pattern, Substitution

    w_ = WildDot("w_")
    w__ = WildPlus("w__")
    w___ = WildStar("w___")

    expr = x + y
    pattern = x + w_
    p, subst = _get_first_match(expr, pattern)
    assert p == Pattern(pattern)
    assert subst == Substitution({'w_': y})

    expr = x + y + z
    pattern = x + w__
    p, subst = _get_first_match(expr, pattern)
    assert p == Pattern(pattern)
    assert subst == Substitution({'w__': Multiset([y, z])})

    expr = x + y + z
    pattern = x + y + z + w___
    p, subst = _get_first_match(expr, pattern)
    assert p == Pattern(pattern)
    assert subst == Substitution({'w___': Multiset()})
