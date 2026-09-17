from sympy import Ne
from sympy.abc import x
from sympy.core.singleton import S
from sympy.solvers.inequalities import reduce_inequalities


def test_reduce_inequalities_preserves_singularities():
    assert reduce_inequalities(1/x <= 1/x, x) == Ne(x, 0)

    e = x/(x - 1) + 1/x <= x + 1/x
    rv = reduce_inequalities(e, x)
    assert rv.subs(x, 0) is S.false
    assert rv.subs(x, S.Half) is S.true
    assert rv.subs(x, 1) is S.false
    assert rv.subs(x, 2) is S.true
