from __future__ import annotations

from sympy.assumptions.ask import Q
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.reasoning_engine import ReasoningEngine
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.testing.pytest import raises


x, y, z = symbols('x y z', real=True)


def test_reasoning_engine_lra():
    for proposition, assumptions, expected in [
        (Q.gt(x, z), Q.gt(x, y) & Q.gt(y, z), True),
        (Q.le(x, z), Q.gt(x, y) & Q.gt(y, z), False),
        (Q.gt(x, y), S.true, None),
        (Q.gt(2, 1), S.true, True),
        (Q.lt(2, 1), S.true, False),
        (S.true, Q.gt(x, 0), True),
        (S.false, Q.gt(x, 0), False),
    ]:
        factbase = EncodedCNF()
        factbase.add_prop(assumptions)
        engine = ReasoningEngine(factbase, use_lra_theory=True)
        query = engine.create_query(CNF.from_prop(proposition),
                                    CNF.from_prop(~proposition))
        assert engine.ask_query(query) is expected


def test_reasoning_engine_lra_inconsistent():
    factbase = EncodedCNF()
    factbase.add_prop(Q.gt(x, y) & Q.gt(y, z) & Q.gt(z, x))
    engine = ReasoningEngine(factbase, use_lra_theory=True)
    with raises(ValueError, match='Inconsistent assumptions'):
        query = engine.create_query(CNF.from_prop(Q.gt(x, 0)),
                                    CNF.from_prop(Q.le(x, 0)))
        engine.ask_query(query)
