from sympy.logic.algorithms.forward_chaining_theory import FCSolver
from sympy.assumptions import Q
from sympy.abc import x, y
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.logic.boolalg import And


fc = FCSolver()
def test_fc_theory():
    fc.reset_state()
    preds = [~Q.positive(x)]
    res = fc.check(preds)
    assert res[0] is True

    fc.reset_state()
    preds = [~Q.positive(x), Q.prime(x)]
    res = fc.check(preds)
    assert res[0] is False
    assert res[1] == [Q.positive(x), ~Q.prime(x)]

    fc.reset_state()
    res = fc.check(preds)
    preds = [~Q.positive(x), Q.prime(x), Q.finite(x)]
    assert res[0] is False
    assert res[1] == [Q.positive(x), ~Q.prime(x)]

    fc.reset_state()
    preds = [~Q.positive(x)]
    res = fc.check(preds)
    assert res[0] is True


def test_multiple_variables():
    fc.reset_state()
    preds = [Q.prime(y), ~Q.negative(x), ~Q.prime(x)]
    res = fc.check(preds)
    assert res[0] is True

    fc.reset_state()
    preds = [Q.positive(x), Q.negative(y)]
    res = fc.check(preds)
    assert res[0] is True

    fc.reset_state()
    preds = [~Q.negative(x), ~Q.prime(x)]
    res = fc.check(preds)
    assert res[0] is True

    fc.reset_state()
    preds = [Q.prime(y), ~Q.negative(x), ~Q.prime(x)]
    res = fc.check(preds)
    assert res[0] is True





