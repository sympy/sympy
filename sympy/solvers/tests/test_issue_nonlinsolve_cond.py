from sympy import symbols, log, Eq, ConditionSet, S
from sympy.solvers.solveset import nonlinsolve

def test_nonlinsolve_conditionset_handling():
    x, y = symbols('x y', real=True)
    # Regression test for incorrect handling of ConditionSet in nonlinsolve
    system = [-y + 4*log(x), -4 + 4*log(y)/x]
    res = nonlinsolve(system, [x, y])

    cond_eq = Eq(-y + 4*log(log(y)), 0)
    cond_set_complex = ConditionSet(y, cond_eq, S.Complexes)
    cond_set_real = ConditionSet(y, cond_eq, S.Reals)

    expected_elem1 = (log(y), cond_set_complex)
    expected_elem2 = (log(y), cond_set_real)

    # We check if these are in the result.
    # Result might contain other permutations or duplicate/similar items,
    # but at least these should be present.
    # Note: res is a FiniteSet
    assert expected_elem1 in res
    assert expected_elem2 in res
