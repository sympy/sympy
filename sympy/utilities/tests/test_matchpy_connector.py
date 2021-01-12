from sympy import symbols, cos, sin
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


def test_matchpy_optional():
    if matchpy is None:
        return

    from matchpy import Pattern, Substitution
    from matchpy import ManyToOneReplacer, ReplacementRule

    p = WildDot("p", optional=1)
    q = WildDot("q", optional=0)

    pattern = p*x + q

    expr1 = 2*x
    pa, subst = _get_first_match(expr1, pattern)
    assert pa == Pattern(pattern)
    assert subst == Substitution({'p': 2, 'q': 0})

    expr2 = x + 3
    pa, subst = _get_first_match(expr2, pattern)
    assert pa == Pattern(pattern)
    assert subst == Substitution({'p': 1, 'q': 3})

    expr3 = x
    pa, subst = _get_first_match(expr3, pattern)
    assert pa == Pattern(pattern)
    assert subst == Substitution({'p': 1, 'q': 0})

    expr4 = x*y + z
    pa, subst = _get_first_match(expr4, pattern)
    assert pa == Pattern(pattern)
    assert subst == Substitution({'p': y, 'q': z})

    replacer = ManyToOneReplacer()
    replacer.add(ReplacementRule(Pattern(pattern), lambda p, q: sin(p)*cos(q)))
    assert replacer.replace(expr1) == sin(2)*cos(0)
    assert replacer.replace(expr2) == sin(1)*cos(3)
    assert replacer.replace(expr3) == sin(1)*cos(0)
    assert replacer.replace(expr4) == sin(y)*cos(z)
