from sympy.core.expr import unchanged
from sympy.core.singleton import S
from sympy.sets.powerset import PowerSet
from sympy.sets.sets import FiniteSet
from sympy.utilities.pytest import raises, XFAIL


def test_powerset_creation():
    assert PowerSet(FiniteSet(1, 2)) == \
        FiniteSet(S.EmptySet, FiniteSet(1), FiniteSet(2), FiniteSet(1, 2))
    assert PowerSet(S.EmptySet) == FiniteSet(S.EmptySet)
    raises(ValueError, lambda: PowerSet(123))
    assert unchanged(PowerSet, S.Reals)
    assert unchanged(PowerSet, S.Integers)


def test_powerset__contains__():
    subset_series = [
        S.EmptySet,
        FiniteSet(1, 2),
        S.Naturals,
        S.Naturals0,
        S.Integers,
        S.Rationals,
        S.Reals,
        S.Complexes]

    l = len(subset_series)
    for i in range(l):
        for j in range(l):
            try:
                if i <= j:
                    assert subset_series[i] in PowerSet(subset_series[j])
                else:
                    assert subset_series[i] not in PowerSet(subset_series[j])
            except:
                raise AssertionError(
                    'Powerset membership test failed between '
                    '{} and {}.'.format(subset_series[i], subset_series[j]))


def test_powerset_contains():
    subset_series = [
        S.EmptySet,
        FiniteSet(1, 2),
        S.Naturals,
        S.Naturals0,
        S.Integers,
        S.Rationals,
        S.Reals,
        S.Complexes]

    l = len(subset_series)
    for i in range(l):
        for j in range(l):
            if i <= j:
                subset_series[i] in PowerSet(subset_series[j])