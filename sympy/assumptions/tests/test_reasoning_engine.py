from __future__ import annotations

from sympy.assumptions.ask import Q
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.reasoning_engine import ReasoningEngine
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.testing.pytest import raises


x, y, z = symbols('x y z', real=True)


def _ask_lra(proposition, assumptions=S.true):
    factbase = EncodedCNF()
    factbase.add_prop(assumptions)
    engine = ReasoningEngine(factbase, use_lra_theory=True)
    query = engine.create_query(CNF.from_prop(proposition),
                                CNF.from_prop(~proposition))
    return engine.ask_query(query)


def test_reasoning_engine_lra():
    assumptions = Q.gt(x, y) & Q.gt(y, z)
    assert _ask_lra(Q.gt(x, z), assumptions) is True
    assert _ask_lra(Q.le(x, z), assumptions) is False
    assert _ask_lra(Q.gt(x, y)) is None
    assert _ask_lra(Q.gt(2, 1)) is True
    assert _ask_lra(Q.lt(2, 1)) is False
    assert _ask_lra(S.true, Q.gt(x, 0)) is True
    assert _ask_lra(S.false, Q.gt(x, 0)) is False


def test_reasoning_engine_lra_inconsistent():
    assumptions = Q.gt(x, y) & Q.gt(y, z) & Q.gt(z, x)
    with raises(ValueError, match='Inconsistent assumptions'):
        _ask_lra(Q.gt(x, 0), assumptions)
