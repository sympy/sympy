from __future__ import annotations

from sympy.assumptions.ask import Q
from sympy.assumptions.cnf import CNF
from sympy.assumptions.lra_satask import check_satisfiability, lra_satask
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.logic.algorithms.lra_theory import UnhandledInput
from sympy.testing.pytest import raises


x, y = symbols('x y', real=True)


def test_cnf_inputs_are_not_mutated():
    positive = symbols('positive', positive=True)
    prop = CNF.from_prop(Q.gt(positive, 0))
    negated = CNF.from_prop(~Q.gt(positive, 0))
    facts = CNF.from_prop(Q.eq(x, 0))
    originals = [cnf.copy() for cnf in (prop, negated, facts)]
    assert check_satisfiability(prop, negated, facts) is True
    assert [cnf.clauses for cnf in (prop, negated, facts)] == [
        cnf.clauses for cnf in originals]


def test_cnf_disjunction_rewriting():
    facts = Q.eq(x, 0) & Q.eq(y, 0)
    for pred in (Q.ne(x, 0), ~Q.eq(x, 0), Q.nonzero(x)):
        assert lra_satask(pred | Q.gt(y, 0), facts) is False
        assert lra_satask(~pred & Q.le(y, 0), facts) is True
    assert lra_satask(S.true, facts) is True
    assert lra_satask(S.false, facts) is False
    with raises(ValueError, match='Inconsistent assumptions'):
        lra_satask(Q.gt(x, 0), S.false)


def test_cnf_known_real_predicates():
    a = symbols('a')
    facts = CNF.from_prop(Q.real(a) & Q.gt(a, 1))
    for pred in (Q.real(a), ~Q.positive_infinite(a), Q.positive(a)):
        prop, negated = CNF.from_prop(pred), CNF.from_prop(~pred)
        assert check_satisfiability(prop, negated, facts, {a}) is True
        with raises(UnhandledInput):
            check_satisfiability(prop, negated, facts)
    with raises(UnhandledInput):
        lra_satask(Q.gt(x, 0), Q.extended_nonnegative(x))
    with raises(UnhandledInput):
        lra_satask(Q.gt(x, 0), Q.prime(x))
