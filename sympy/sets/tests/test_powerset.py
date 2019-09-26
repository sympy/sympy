from sympy.core.expr import unchanged
from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.sets.contains import Contains
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
                    assert subset_series[i] in \
                        PowerSet(subset_series[j], evaluate=False)
                else:
                    assert subset_series[i] not in \
                        PowerSet(subset_series[j], evaluate=False)
            except:
                raise AssertionError(
                    'Powerset membership test failed between '
                    '{} and {}.'.format(subset_series[i], subset_series[j]))


@XFAIL
def test_failing_powerset__contains__():
    # XXX These are failing when evaluate=True,
    # but using unevaluated PowerSet works fine.
    assert FiniteSet(1, 2) not in PowerSet(S.EmptySet)
    assert S.Naturals not in PowerSet(S.EmptySet)
    assert S.Naturals not in PowerSet(FiniteSet(1, 2))
    assert S.Naturals0 not in PowerSet(S.EmptySet)
    assert S.Naturals0 not in PowerSet(FiniteSet(1, 2))
    assert S.Integers not in PowerSet(S.EmptySet)
    assert S.Integers not in PowerSet(FiniteSet(1, 2))
    assert S.Rationals not in PowerSet(S.EmptySet)
    assert S.Rationals not in PowerSet(FiniteSet(1, 2))
    assert S.Reals not in PowerSet(S.EmptySet)
    assert S.Reals not in PowerSet(FiniteSet(1, 2))
    assert S.Complexes not in PowerSet(S.EmptySet)
    assert S.Complexes not in PowerSet(FiniteSet(1, 2))


def test_powerset__len__():
    A = PowerSet(S.EmptySet, evaluate=False)
    assert len(A) == 1
    A = PowerSet(A, evaluate=False)
    assert len(A) == 2
    A = PowerSet(A, evaluate=False)
    assert len(A) == 4
    A = PowerSet(A, evaluate=False)
    assert len(A) == 16


def test_powerset_contains():
    A = PowerSet(FiniteSet(1), evaluate=False)
    assert A.contains(2) == Contains(2, A)

    x = Symbol('x')

    A = PowerSet(FiniteSet(x), evaluate=False)
    assert A.contains(FiniteSet(1)) == Contains(A, FiniteSet(1))
