from sympy.abc import x
from sympy.sets.fancysets import Interval
from sympy.sets.sets import Set, FiniteSet, Union
from sympy.sets.indicator import indicator


def test_indicator_function():
    A = Interval(0, 1)
    assert indicator(0.5, A) == 1
    assert indicator(-0.5, A) == 0
    B = Union(Interval(0, 1), Interval(2, 3))
    assert indicator(0.5, B) == 1
    assert indicator(1.5, B) == 0
    assert indicator(2.5, B) == 1
    C = FiniteSet(1, 2, 3)
    assert indicator(0, C) == 0
    assert indicator(1, C) == 1


def test_indicator_symbolic_function():
    A = Set('A')
    interval = Interval(0, 1)
    symbolic_indicator = indicator(x, A)
    assert symbolic_indicator.subs({x: 1}) == indicator(1, A)
    assert symbolic_indicator.subs({A: interval}) == indicator(x, interval)
