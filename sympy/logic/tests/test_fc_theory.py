from sympy.logic.algorithms.forward_chaining_theory import FCSolver
from sympy.assumptions import Q
from sympy.abc import x, y
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.logic.boolalg import And

data = []
encoding = {
    Q.positive(x): 1,
    Q.prime(x): 2,
    Q.finite(x): 3,
    Q.odd(x) : 4,
    Q.even(x): 5,
    Q.integer(x): 6,
    Q.commutative(x): 7,
}
cnf = EncodedCNF(data, encoding)
fc, state = FCSolver.from_encoded_cnf(cnf, testing_mode=True)

def test_check():

    fc.reset_state()
    preds = [6, -5, -4]
    res = fc.check(preds, False)
    assert res[0] is False

    fc.reset_state()
    preds = [-1]
    res = fc.check(preds, False)
    assert res[0] is True

    fc.reset_state()
    preds = [-1, 2]
    res = fc.check(preds, False)
    assert res[0] is False
    assert res[1] == [-2, 1]

    fc.reset_state()
    preds = [-1, 2, 3]
    res = fc.check(preds, False)
    assert res[0] is False
    assert res[1] == [-2, 1]

    fc.reset_state()
    preds = [-1]
    res = fc.check(preds, False)
    assert res[0] is True

    fc.reset_state()
    preds = [6, -5, -4]
    res = fc.check(preds, False)
    assert res[0] is False

    # satask(Q.odd(x * y), Q.odd(x) & Q.odd(y)) is True


def test_immediate_conflict_detection():

    assert not fc.lit_in_theory(0)
    assert fc.lit_in_theory(1)
    assert fc.lit_in_theory(2)
    assert fc.lit_in_theory(3)
    assert fc.lit_in_theory(-1)
    assert fc.lit_in_theory(-2)
    assert fc.lit_in_theory(-3)

    fc.reset_state()
    lits = [-1]
    res = (True, None)
    for lit in lits:
        res = fc.assert_lit(lit, state)
        if not res[0]:
            break
    assert res[0] is True

    fc.reset_state()
    lits = [-1, 2]
    res = (True, None)
    for lit in lits:
        res = fc.assert_lit(lit, state)
        if not res[0]:
            break
    assert res[0] is False

#
# def test_multiple_variables():
#     fc.reset_state()
#     preds = [Q.prime(y), ~Q.negative(x), ~Q.prime(x)]
#     res = fc.check(preds)
#     assert res[0] is True
#
#     fc.reset_state()
#     preds = [Q.positive(x), Q.negative(y)]
#     res = fc.check(preds)
#     assert res[0] is True
#
#     fc.reset_state()
#     preds = [~Q.negative(x), ~Q.prime(x)]
#     res = fc.check(preds)
#     assert res[0] is True
#
#     fc.reset_state()
#     preds = [Q.prime(y), ~Q.negative(x), ~Q.prime(x)]
#     res = fc.check(preds)
#     assert res[0] is True





